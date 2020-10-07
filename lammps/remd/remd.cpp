/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Trung Nguyen (Northwestern University)
   Reference: Bunker and Dunweg, Phys. Rev. E, 63, 016701, 2000. 
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "remd.h"
#include "universe.h"
#include "domain.h"
#include "atom.h"
#include "update.h"
#include "integrate.h"
#include "modify.h"
#include "compute.h"
#include "force.h"
#include "kspace.h"
#include "input.h"
#include "variable.h"
#include "output.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "thermo.h"
#include "fix.h"
#include "random_park.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

// adopt from fix_adapt for parsing pair/atom/kspace from arguments

//#define REMD_DEBUG 1

enum{PAIR,DIELECTRIC,KSPACE,ATOM};

/* ---------------------------------------------------------------------- */

REMD::REMD(LAMMPS *lmp) : Pointers(lmp), nexchange(0), exchange(NULL)
{
  fixgpu = NULL;
}

/* ---------------------------------------------------------------------- */

REMD::~REMD()
{
  for (int m = 0; m < nexchange; m++) {
    if (exchange[m].which == PAIR) {
      delete [] exchange[m].pstyle;
      delete [] exchange[m].pparam;
    }
  }

  delete [] exchange;

  MPI_Comm_free(&roots);
  if (ranswap) delete ranswap;
  delete ranboltz;
  delete [] set_param;
  delete [] param2world;
  delete [] world2param;
  delete [] world2root;
}

/* ----------------------------------------------------------------------
   perform tempering with inter-world swaps
------------------------------------------------------------------------- */

void REMD::command(int narg, char **arg)
{
  if (universe->nworlds == 1)
    error->all(FLERR,"Must have more than one processor partition to remd");
  if (domain->box_exist == 0)
    error->all(FLERR,"remd command before simulation box is defined");

  int nsteps = force->inumeric(FLERR,arg[0]);
  nevery = force->inumeric(FLERR,arg[1]);
  temperature = force->numeric(FLERR,arg[2]);
  param = force->numeric(FLERR,arg[3]);
  seed_swap = force->inumeric(FLERR,arg[4]);
  seed_boltz = force->inumeric(FLERR,arg[5]);

  my_set_param = universe->iworld;
  exchange = 0;

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"pair") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal remd command");
      nexchange++;
      iarg += 5;
    } else if (strcmp(arg[iarg],"dielectric") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal remd command");
      nexchange++;
      iarg += 1;
    } else break;
  }

  if (nexchange == 0) error->all(FLERR,"Illegal remd command");
  exchange = new Exchange[nexchange];

  iarg = 6;
  nexchange = 0;
  int restart = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"pair") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal remd command");
      exchange[nexchange].which = PAIR;
      int n = strlen(arg[iarg+1]) + 1;
      exchange[nexchange].pstyle = new char[n];
      strcpy(exchange[nexchange].pstyle,arg[iarg+1]);
      n = strlen(arg[iarg+2]) + 1;
      exchange[nexchange].pparam = new char[n];
      exchange[nexchange].pair = NULL;
      strcpy(exchange[nexchange].pparam,arg[iarg+2]);
      force->bounds(FLERR,arg[iarg+3],atom->ntypes,
                    exchange[nexchange].ilo,exchange[nexchange].ihi);
      force->bounds(FLERR,arg[iarg+4],atom->ntypes,
                    exchange[nexchange].jlo,exchange[nexchange].jhi);
      
      nexchange++;
      iarg += 5;
    } else if (strcmp(arg[iarg],"dielectric") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal remd command");
      exchange[nexchange].which = DIELECTRIC;

      if (force->kspace != NULL) {
        kspace_qqrd2e = (double *) force->kspace->extract("qqrd2e");
      }
      
      nexchange++;
      iarg += 1;
    } else if (strcmp(arg[iarg],"restart") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal remd command");
      my_set_param = force->inumeric(FLERR,arg[iarg+1]);
      if ((my_set_param < 0) || (my_set_param >= universe->nworlds))
        error->universe_one(FLERR,"Illegal param index");
      restart = 1;
      iarg += 2;
    } else break;
  }

  // swap frequency must evenly divide total # of timesteps

  if (nevery <= 0)
    error->universe_all(FLERR,"Invalid frequency in remd command");
  nswaps = nsteps/nevery;
  if (nswaps*nevery != nsteps)
    error->universe_all(FLERR,"Non integer # of swaps in remd command");

  // setup for long tempering run

  update->whichflag = 1;
  update->nsteps = nsteps;
  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->firststep + nsteps;
  if (update->laststep < 0)
    error->all(FLERR,"Too many timesteps");

  lmp->init();

  // local storage

  me_universe = universe->me;
  MPI_Comm_rank(world,&me);
  nworlds = universe->nworlds;
  iworld = universe->iworld;
  boltz = force->boltz;

  // create MPI communicator for root proc from each world

  int color;
  if (me == 0) color = 0;
  else color = 1;
  MPI_Comm_split(universe->uworld,color,0,&roots);

  // RNGs for swaps and Boltzmann test
  // warm up Boltzmann RNG

  if (seed_swap) ranswap = new RanPark(lmp,seed_swap);
  else ranswap = NULL;
  ranboltz = new RanPark(lmp,seed_boltz + me_universe);
  for (int i = 0; i < 100; i++) ranboltz->uniform();

  // world2root[i] = global proc that is root proc of world i

  world2root = new int[nworlds];
  if (me == 0)
    MPI_Allgather(&me_universe,1,MPI_INT,world2root,1,MPI_INT,roots);
  MPI_Bcast(world2root,nworlds,MPI_INT,0,world);

  // create static list of set params
  // allgather tempering arg "temp" across root procs
  // bcast from each root to other procs in world

  set_param = new double[nworlds];
  if (me == 0) MPI_Allgather(&param,1,MPI_DOUBLE,set_param,1,MPI_DOUBLE,roots);
  MPI_Bcast(set_param,nworlds,MPI_DOUBLE,0,world);

  // create world2temp only on root procs from my_set_temp
  // create temp2world on root procs from world2temp,
  //   then bcast to all procs within world

  world2param = new int[nworlds];
  param2world = new int[nworlds];
  if (me == 0) {
    MPI_Allgather(&my_set_param,1,MPI_INT,world2param,1,MPI_INT,roots);
    for (int i = 0; i < nworlds; i++) param2world[world2param[i]] = i;
  }
  MPI_Bcast(param2world,nworlds,MPI_INT,0,world);

  // if restarting remd, reset parameter

  if (restart) {
    exchange_settings(set_param[my_set_param]);
  }

  // setup tempering runs

  int i,which,partner,swap,partner_set_param,partner_world,eflag,vflag;
  double pe[2],pe_partner[2],boltz_factor,new_param;

  if (me_universe == 0 && universe->uscreen)
    fprintf(universe->uscreen,"Setting up tempering ...\n");

  update->integrate->setup(1);

  if (me_universe == 0) {
    if (universe->uscreen) {
      fprintf(universe->uscreen,"Step");
      for (int i = 0; i < nworlds; i++)
        fprintf(universe->uscreen," T%d",i);
      fprintf(universe->uscreen,"\n");
    }
    if (universe->ulogfile) {
      fprintf(universe->ulogfile,"Step");
      for (int i = 0; i < nworlds; i++)
        fprintf(universe->ulogfile," T%d",i);
      fprintf(universe->ulogfile,"\n");
    }
    print_status();
  }

  timer->init();
  timer->barrier_start();

  // initialize the exchanged parameter(s)

  init();
  eflag = 1;
  vflag = 0;

  for (int iswap = 0; iswap < nswaps; iswap++) {

    // run for nevery timesteps

    update->integrate->run(nevery);

    // compute PE

    if (force->pair && force->pair->compute_flag) {
      force->pair->compute(eflag,vflag);
    }

    // accumulate force/energy/virial from /gpu pair styles
    if (fixgpu) fixgpu->post_force(vflag);

    // pe[0] -> H1(x)

    pe[0] = compute_epair();

    // which = which of 2 kinds of swaps to do (0,1)

    if (!ranswap) which = iswap % 2;
    else if (ranswap->uniform() < 0.5) which = 0;
    else which = 1;

    // partner_set_temp = which set temp I am partnering with for this swap

    if (which == 0) {
      if (my_set_param % 2 == 0) partner_set_param = my_set_param + 1;
      else partner_set_param = my_set_param - 1;
    } else {
      if (my_set_param % 2 == 1) partner_set_param = my_set_param + 1;
      else partner_set_param = my_set_param - 1;
    }

    // partner = proc ID to swap with
    // if partner = -1, then I am not a proc that swaps

    partner = -1;
    if (partner_set_param >= 0 && partner_set_param < nworlds) {
      partner_world = param2world[partner_set_param];
      partner = world2root[partner_world];
    }

    // if this world is to be exchanged,
    // compute potential energy with the parameter from partner 

    if (partner != -1) {

      exchange_settings(set_param[partner_set_param]);

      if (force->pair && force->pair->compute_flag) {
        force->pair->compute(eflag,vflag);
      }

      // accumulate force/energy/virial from /gpu pair styles
      if (fixgpu) fixgpu->post_force(vflag);

      // pe[1] -> H2(x)
      pe[1] = compute_epair();
    }

    // swap with a partner, only root procs in each world participate
    // hi proc sends PE to low proc
    // lo proc make Boltzmann decision on whether to swap
    // lo proc communicates decision back to hi proc

    swap = 0;
    if (me == 0 && partner != -1) {
      if (me_universe > partner)
        MPI_Send(&pe,2,MPI_DOUBLE,partner,0,universe->uworld);
      else
        MPI_Recv(&pe_partner,2,MPI_DOUBLE,partner,0,universe->uworld,MPI_STATUS_IGNORE);

      // boltz_factor from Eq. (6) in Ref.
      // pe[0] -> H1(x)
      // pe[1] -> H2(x)
      // pe_partner[0] -> H2(y)
      // pe_partner[1] -> H1(y)

      if (me_universe < partner) {
        double beta = 1.0/(boltz*temperature);
        boltz_factor = -beta*(pe[1] + pe_partner[1] - pe[0] - pe_partner[0]);
        if (boltz_factor >= 0.0) swap = 1;
        else if (ranboltz->uniform() < exp(boltz_factor)) swap = 1;

      }

      if (me_universe < partner)
        MPI_Send(&swap,1,MPI_INT,partner,0,universe->uworld);
      else
        MPI_Recv(&swap,1,MPI_INT,partner,0,universe->uworld,MPI_STATUS_IGNORE);

#ifdef REMD_DEBUG
      if (me_universe < partner)
        printf("SWAP %d & %d: yes = %d, params = %d %d, H1(x) = %g H2(x) = %g; "
               " H2(y) = %g H1(y) = %g, Bz = %g %g\n", me_universe,partner,swap,
               my_set_param,partner_set_param, pe[0],pe[1],pe_partner[0],pe_partner[1],
               boltz_factor,exp(boltz_factor));
#endif

    }

    // bcast swap result to other procs in my world

    MPI_Bcast(&swap,1,MPI_INT,0,world);


    // if my world is not swapped, restore the parameter

    if (!swap) {

      exchange_settings(set_param[my_set_param]);
    }

    // update my_set_temp and temp2world on every proc
    // root procs update their value if swap took place
    // allgather across root procs
    // bcast within my world

    if (swap) my_set_param = partner_set_param;
    if (me == 0) {
      MPI_Allgather(&my_set_param,1,MPI_INT,world2param,1,MPI_INT,roots);
      for (i = 0; i < nworlds; i++) param2world[world2param[i]] = i;
    }
    MPI_Bcast(param2world,nworlds,MPI_INT,0,world);

    // print out current swap status

    if (me_universe == 0) print_status();
  }

  timer->barrier_stop();

  update->integrate->cleanup();

  Finish finish(lmp);
  finish.end(1);

  update->whichflag = 0;
  update->firststep = update->laststep = 0;
  update->beginstep = update->endstep = 0;
}

/* ---------------------------------------------------------------------- */

void REMD::init()
{
  int i,j;

  // setup and error checks

  anypair = 0;

  for (int m = 0; m < nexchange; m++) {
    Exchange *exch = &exchange[m];

    if (exch->which == PAIR) {
      anypair = 1;
      exch->pair = NULL;

      // if ad->pstyle has trailing sub-style annotation ":N",
      //   strip it for pstyle arg to pair_match() and set nsub = N
      // this should work for appended suffixes as well

      int n = strlen(exch->pstyle) + 1;
      char *pstyle = new char[n];
      strcpy(pstyle,exch->pstyle);

      char *cptr;
      int nsub = 0;
      if ((cptr = strchr(pstyle,':'))) {
        *cptr = '\0';
        nsub = force->inumeric(FLERR,cptr+1);
      }

      if (lmp->suffix_enable) {
        int len = 2 + strlen(pstyle) + strlen(lmp->suffix);
        char *psuffix = new char[len];
        strcpy(psuffix,pstyle);
        strcat(psuffix,"/");
        strcat(psuffix,lmp->suffix);
        exch->pair = force->pair_match(psuffix,1,nsub);
        delete[] psuffix;
      }
      if (exch->pair == NULL) exch->pair = force->pair_match(pstyle,1,nsub);
      if (exch->pair == NULL)
        error->all(FLERR,"remd pair style does not exist");

      void *ptr = exch->pair->extract(exch->pparam,exch->pdim);
      if (ptr == NULL)
        error->all(FLERR,"remd pair style param not supported");

      // for pair styles only parameters that are 2-d arrays in atom types or
      // scalars are supported

      if (exch->pdim != 2 && exch->pdim != 0)
        error->all(FLERR,"remd pair style param is not compatible");

      if (exch->pdim == 2) exch->array = (double **) ptr;
      if (exch->pdim == 0) exch->scalar = (double *) ptr;

      // if pair hybrid, test that ilo,ihi,jlo,jhi are valid for sub-style

      if (strcmp(force->pair_style,"hybrid") == 0 ||
          strcmp(force->pair_style,"hybrid/overlay") == 0) {
        PairHybrid *pair = (PairHybrid *) force->pair;
        for (i = exch->ilo; i <= exch->ihi; i++)
          for (j = MAX(exch->jlo,i); j <= exch->jhi; j++)
            if (!pair->check_ijtype(i,j,pstyle))
              error->all(FLERR,"remd type pair range is not valid for "
                         "pair hybrid sub-style");
      }

      delete [] pstyle;
    }
  }

  // detect if package gpu is present

  int ifixgpu = modify->find_fix("package_gpu");
  if (ifixgpu >= 0) fixgpu = modify->fix[ifixgpu];
}

/* ----------------------------------------------------------------------
   exchange pair parameters
------------------------------------------------------------------------- */

void REMD::exchange_settings(double value)
{
  int i,j;

  for (int m = 0; m < nexchange; m++) {
    Exchange *exch = &exchange[m];

    // set global scalar or type pair array values

    if (exch->which == PAIR) {
      if (exch->pdim == 0) {
        *exch->scalar = value;
      } else if (exch->pdim == 2) {
        for (i = exch->ilo; i <= exch->ihi; i++)
          for (j = MAX(exch->jlo,i); j <= exch->jhi; j++)
            exch->array[i][j] = value;
      }
    } else if (exch->which == DIELECTRIC) {

      force->dielectric = value;
      force->qqrd2e = force->qqr2e/value;  // pair
      if (force->kspace) *kspace_qqrd2e = force->qqrd2e;      // kspace
    }
  }

  // re-initialize pair styles if any PAIR settings were changed
  // this resets other coeffs that may depend on changed values,
  // and also offset and tail corrections
  if (anypair) {
    for (int m = 0; m < nexchange; m++) {
      Exchange *exch = &exchange[m];
      if (exch->which == PAIR) {
        exch->pair->reinit();
      }
    }
  }

}

/* ----------------------------------------------------------------------
   proc 0 prints current tempering status
------------------------------------------------------------------------- */

void REMD::print_status()
{
  if (universe->uscreen) {
    fprintf(universe->uscreen,BIGINT_FORMAT,update->ntimestep);
    for (int i = 0; i < nworlds; i++)
      fprintf(universe->uscreen," %d",world2param[i]);
    fprintf(universe->uscreen,"\n");
  }
  if (universe->ulogfile) {
    fprintf(universe->ulogfile,BIGINT_FORMAT,update->ntimestep);
    for (int i = 0; i < nworlds; i++)
      fprintf(universe->ulogfile," %d",world2param[i]);
    fprintf(universe->ulogfile,"\n");
    fflush(universe->ulogfile);
  }
}

/* ----------------------------------------------------------------------
   obtain pair energy from lammps accumulators
------------------------------------------------------------------------- */

double REMD::compute_epair()
{
  double eng, eng_pair;
  eng = 0.0;
  if (force->pair)
    eng = force->pair->eng_vdwl + force->pair->eng_coul;
  MPI_Allreduce(&eng,&eng_pair,1,MPI_DOUBLE,MPI_SUM,world);

  return eng_pair;
}

