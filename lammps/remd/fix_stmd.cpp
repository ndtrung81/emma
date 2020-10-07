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

#include "string.h"
#include "fix_stmd.h"
#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "update.h"
#include "memory.h"
#include "modify.h"
#include "error.h"
#include "force.h"
#include <math.h>

using namespace LAMMPS_NS;
using namespace FixConst;

#define EPSILON 1.0e-2

/* ----------------------------------------------------------------------
   Contributing author: Trung Dac Nguyen (ORNL/Northwestern)
------------------------------------------------------------------------- */

FixSTMD::FixSTMD(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 10) error->all(FLERR,"Illegal fix stmd command");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 5;
  extscalar = 0;
  extvector = 0;
  restart_global = 1;
  global_freq = 1;

  f_current = 1.0+2.5e-4;
  f_final = 1.0+1.0e-8;
  flatcheck = 200000;
  nbins_left_edge = 2;
  nbins_right_edge = 2;
  threshold = 0.2;
  fraction = fraction0 = 0.6;
  maxcount = 1000000;
  smooth = 2;
  padflag = 0;
  fixed = 1;
  ke_included = 0;

  char* defaultname = (char*)"hist_*.txt";
  int n = strlen(defaultname) + 1;
  filename = new char[n];
  strcpy(filename,defaultname);
  
  int tflag = 0;
  int eflag = 0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"temp") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix stmd command");
      tflag = 1;
      T_l = force->numeric(FLERR,arg[iarg+1]);
      T_h = force->numeric(FLERR,arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"energy") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix stmd command");
      eflag = 1;
      E_min = force->numeric(FLERR,arg[iarg+1]);
      E_max = force->numeric(FLERR,arg[iarg+2]);
      binsize_E = force->numeric(FLERR,arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"fd") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix stmd command");
      double fd = force->numeric(FLERR,arg[iarg+1]);
      if (fd < 0.0) error->all(FLERR,"fix stmd: initial fd should be positive");
      f_current = 1.0 + fd;
      fd = force->numeric(FLERR,arg[iarg+2]);
      if (fd < 0.0) error->all(FLERR,"fix stmd: final fd should be positive");
      f_final = 1.0 + fd;
      iarg += 3;
    } else if (strcmp(arg[iarg],"flatcheck") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix stmd command");
      flatcheck = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"left_edge") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix stmd command");
      nbins_left_edge = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"right_edge") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix stmd command");
      nbins_right_edge = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"threshold") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix stmd command");
      threshold = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"smooth") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix stmd command");
      smooth = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"fraction") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix stmd command");
      fraction0 = force->numeric(FLERR,arg[iarg+1]);
      fraction = fraction0;
      iarg += 2;
    } else if (strcmp(arg[iarg],"maxcount") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix stmd command");
      maxcount = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix stmd command");
      delete [] filename;
      n = strlen(arg[iarg+1]) + 1;
      filename = new char[n];
      strcpy(filename,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"pad") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix stmd command");
      padflag = force->inumeric(FLERR,arg[iarg+1]);
      if (padflag < 0) error->all(FLERR,"Illegal fix stmd command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"fixed") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix stmd command");
      fixed = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"kinetic") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix stmd command");
      ke_included = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix stmd command");
  }

  // sanity checks

  if (!tflag || !eflag)
    error->all(FLERR,"Illegal fix stmd command without temp or "
                     "energy windows");
  if (T_h < T_l)
    error->all(FLERR,"Illegal fix stmd command: T_h should be "
                     "greater than T_l");
  if (E_max < E_min)
    error->all(FLERR,"Illegal fix stmd command: E_max should be "
                     "greater than E_min");

  if (binsize_E <= 0.0)
    error->all(FLERR,"Illegal fix stmd command: "
                     "energy bin size must be positive");

  // allocate memory for the temperature profile T_E and energy histogram h_E

  nbins = ceil((E_max - E_min) / binsize_E);
  memory->create(T_E,nbins,"fix:t_e");
  memory->create(h_E,nbins,"fix:h_e");
  for (int i = 0; i < nbins; i++) {
    T_E[i] = T_h;
    h_E[i] = 0.0;
  }

  // create the internal energy compute

  n = strlen(id) + 6;
  id_pe = new char[n];
  strcpy(id_pe,id);
  strcat(id_pe,"_pe");

  char **newarg = new char*[3];
  newarg[0] = id_pe;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "pe";
  modify->add_compute(3,newarg);

  if (ke_included) {
    id_temp = new char[n];
    strcpy(id_temp,id);
    strcat(id_temp,"_temp");

    newarg[0] = id_temp;
    newarg[1] = (char *) "all";
    newarg[2] = (char *) "temp";
    modify->add_compute(3,newarg);
  }

  delete [] newarg;

  pe = NULL;
  temperature = NULL;
  flattening = 1;
}

/* ---------------------------------------------------------------------- */

FixSTMD::~FixSTMD()
{
  // delete compute potential energy and kinetic energy

  modify->delete_compute(id_pe);
  delete [] id_pe;
  if (ke_included) {
    modify->delete_compute(id_temp);
    delete [] id_temp;
  }
  
  // deallocate memory
  memory->destroy(T_E);
  memory->destroy(h_E);

  delete [] filename;
}

/* ---------------------------------------------------------------------- */

int FixSTMD::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSTMD::init()
{
  // check if there is a potential energy compute defined

  int icompute = modify->find_compute(id_pe);
  if (icompute < 0)
    error->all(FLERR,"Compute pe for fix stmd does not exist");
  pe = modify->compute[icompute]; 

  if (ke_included) {
    icompute = modify->find_compute(id_temp);
    if (icompute < 0)
      error->all(FLERR,"Compute temp for fix stmd does not exist");
    temperature = modify->compute[icompute]; 
  }

  // check if there is a thermostat is defined

  int thermostatted = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"nvt") == 0 || 
        strcmp(modify->fix[i]->style,"nvt/asphere") == 0 ||
        strcmp(modify->fix[i]->style,"langevin") == 0 ||
        strcmp(modify->fix[i]->style,"temp/berendsen") == 0 ||
        strcmp(modify->fix[i]->style,"temp/csvr") == 0) {
      thermostatted = 1;
      break;
    }

  if (!thermostatted)
    error->warning(FLERR,"fix stmd requires a thermostat to maintain " 
                         "the kinetic temperature at T_h");

  if (f_current < f_final)
    error->warning(FLERR,
      "Current scaling factor is smaller than the final value");
}

/* ---------------------------------------------------------------------- */

void FixSTMD::setup(int vflag)
{
  T_0 = T_h;
  T_current = T_h;
  first_time = 0;
  nout = 0;

  deltaf = force->boltz * log(f_current) / (2.0 * binsize_E);
  post_force(vflag);
}

/* ---------------------------------------------------------------------- 
  Reference for STMD: 
     J. Kim, J. E. Straub, T. Keyes, Phys. Rev. Lett. 97, 050601 (2006).
 ------------------------------------------------------------------------ */

void FixSTMD::post_force(int vflag)
{
  // compute instantaneous system energy

  double E = compute_energy();

  // if E is out of range
  // return if the energy window is fixed and never in
  // otherwise, resize the energy window

  if (E < E_min || E > E_max) {
    if (fixed) {
      nout++;
      if (!first_time) {
        return;
      } else {
        // because the newly generated configuration 
        // is accepted already, I think return makes sense 
        // rather than counting the current E and updating T again,
        // which would lead to a spike often the lower bound
        return;
      }
    } else {
      E_current = E;
      resize_energy_window(E_current);
      first_time = 1;
    }
  } else {
    E_current = E;
    first_time = 1;
  }

  // use E_current for updating histogram h_E and temperature estimate T_E

  int j;
  double alpha, lambdaj, Ej, Tj, gamma, T;

  // find the current bin, neglecting the lower and upper bounds
  // only update h_E (when flattening is done -- previous version)

  j = static_cast<int>((E_current - E_min) / binsize_E);

  if (j == nbins - 1) return;

  h_E[j] += 1.0;

  // update the temperature estimate T_E 
  // see Eq. 3 in Kim et al., PRL 97, 050601, 2006.
  // Important notes from Ref.: 
  // "restricting updates to T_l < T < T_h"
  // "maintaining T = T_l (T_h) beyond the lower 
  // (upper) temperature bounds T_l (T_h)"

  if (j > 0) {
    alpha = 1.0 / (1.0 + deltaf * T_E[j-1]);
    T = alpha * T_E[j-1];
    if (T >= T_l && T <= T_h) T_E[j-1] = T;
    if (j-1 == 0 && flattening == 0) T_E[j-1] = T_l;
  } else {
    T_E[j] = T_l;
  }

  alpha = 1.0 / (1.0 - deltaf * T_E[j+1]);
  T = alpha * T_E[j+1];
  if (T >= T_l && T <= T_h) T_E[j+1] = T;
  if (j+1 == nbins-1) T_E[j+1] = T_h;

  // interpolate for T(E) between [E_j;T_j] and [E_j+1; T_j+1]
  // see Eq. 4
 
  lambdaj = (T_E[j+1] - T_E[j]) / binsize_E;
  Ej = E_min + j * binsize_E;
  Tj = T_E[j] + lambdaj * (E_current - Ej);

  if (Tj > T_h) Tj = T_h;
  else if (Tj < T_l) Tj = T_l;

  T_current = Tj;
  gamma = T_0 / T_current;

  // scale the atom forces

  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      f[i][0] *= gamma;
      f[i][1] *= gamma;
      f[i][2] *= gamma;
    }
  }

  // flattening T_E and flatness checking for h_E

  if (update->ntimestep % flatcheck == 0) {
    if (flattening == 1) low_end_flattening();
    check_flatness_kl();
  }
}

/* ---------------------------------------------------------------------- 
   Flatten T_E for E < E_min when T_min is still greater than T_l
   --------------------------------------------------------------------- */ 

void FixSTMD::low_end_flattening()
{
  int i_min = 0;
  double T_min=T_E[0];
  for (int i = 0; i < nbins; i++) {
    if (T_min > T_E[i]) {
      T_min = T_E[i];
      i_min = i;
    }
  }
    
  if (T_min >= T_l) {
    for (int i = 0; i <= i_min; i++)
      T_E[i] = T_min;
    if (T_min - T_l < EPSILON)
      flattening = 0;
  } 
}

/* ---------------------------------------------------------------------- 
  Check if the energy histogram is flat:
    adjust the scaling factor
    reset the energy histogram 
------------------------------------------------------------------------ */

int FixSTMD::check_flatness()
{
  int i, m, flat;

  // find the edge bins where the jump/drop in H(E) occurs
  // start from the center bin, goes to the left, look for the first drop
  // compare h_E[i] with the running mean for better edge detection

  int center_bin = nbins / 2;
  double left, right;
  double sum, mean, mean_left, mean_right, first_fluctuation, stdev;
  
  left = right = h_E[center_bin];
  left_edge = right_edge = center_bin;

  sum = 0.0;
  m = 0;
  mean_left = left;
  for (i = center_bin-1; i >= 0; i--) {
    if (h_E[i] < (1.0-threshold)*left) {
      left_edge = i+1;
      break;
    } else {
      left_edge = i;
      sum += h_E[i];
      m++;
      mean_left = sum/m;
    }
    left = mean_left; //h_E[i];
  }

  // start from the center bin, goes to the right, look for the first drop
  // compare h_E[i] with the running mean for better edge detection

  sum = 0.0;
  m = 0;
  mean_right = right;
  for (i = center_bin+1; i < nbins; i++) {
    if (h_E[i] < (1.0-threshold)*right) {
      right_edge = i-1;
      break;
    } else {
      right_edge = i;
      sum += h_E[i];
      m++;
      mean_right = sum/m;
    }
    right = mean_right; //h_E[i];
  }

  // to be more conservative by excluding nbins_edge on both sides

  left_edge += nbins_left_edge;
  right_edge -= nbins_right_edge;
  left = h_E[left_edge];
  right = h_E[right_edge];
   
  // check if the energy histogram is flat excluding the edge bins
  // when the number of bins in the checked region is considerable
  // mean = mean value of those between the left and right edges
  // fluctuation = first value greater than the set threshold fraction of the mean

  flat = 0;
  mean = -1.0;
  first_fluctuation = -1.0;
  if ((right_edge - left_edge + 1) >= fraction*nbins && 
      left > 0.0 && right > 0.0) {

    // find the mean value

    m = 0;
    mean = 0.0;
    stdev = 0.0;
    for (i = left_edge; i <= right_edge; i++) {
      mean += h_E[i];
      stdev += h_E[i]*h_E[i];
      m++;
    }
    mean /= (double)m;
    stdev = sqrt(stdev / (double)m - mean * mean);

    flat = 1;
    for (i = left_edge; i <= right_edge; i++) {
      if (fabs((h_E[i] - mean)/mean) > threshold) {
        first_fluctuation = fabs((h_E[i] - mean)/mean);
        flat = 0;
        break;
      }
    }
  } else {
    m = 0;
    mean = 0.0;
    stdev = 0.0;
    for (i = 0; i < nbins; i++) {
      mean += h_E[i];
      stdev += h_E[i]*h_E[i];
      m++;
    }
    mean /= (double)m;
    stdev = sqrt(stdev / (double)m - mean * mean);
  }

  // handle the cases where the fraction value 
  // cannot be reached (fluctuation = -1)
  // after the histogram in the flat region is already sufficiently big
  // find the average value between left and right
  // decrease fraction by 0.05 until it reaches 0.5

  double average = (mean_left + mean_right) * 0.5;
  if (average > maxcount && (right_edge - left_edge + 1) >= 0.5*nbins 
     && flat == 0 && first_fluctuation == -1 && fraction > 0.5) {
    fraction -= 0.05;
  }

  // print out the current histogram h_E and temperature profile T_E

  if (comm->me == 0) {
    write_histograms(mean, first_fluctuation, flat);
  }

  // reduce the scaling factor
  //   if smaller than f_final, then estimate log(g(E)) from T_E
  //   otherwise, reset the energy histogram and fraction

  if (flat == 1) {
    if (f_current <= f_final) {
      dos_estimate();
    } else {
      f_current = sqrt(f_current);
      deltaf = force->boltz * log(f_current) / (2.0 * binsize_E);
      for (i = 0; i < nbins; i++) 
        h_E[i] = 0.0;
      fraction = fraction0;
    }
  }

  return flat;
}

/* ----------------------------------------------------------------------
  Check for histogram flatness based on Kullbach-Leibler
------------------------------------------------------------------------- */
int FixSTMD::check_flatness_kl()
{
  int i, m, flat = 0;

  // find the edge bins where the jump/drop in H(E) occurs
  // start from the center bin, goes to the left, look for the first drop
  // compare h_E[i] with the running mean for better edge detection

  int center_bin = nbins / 2;
  double left, right;
  double sum, mean, mean_left, mean_right, first_fluctuation, stdev;
  
  left = right = h_E[center_bin];
  left_edge = right_edge = center_bin;

  sum = 0.0;
  m = 0;
  mean_left = left;
  for (i = center_bin-1; i >= 0; i--) {
    if (h_E[i] < (1.0-threshold)*left) {
      left_edge = i+1;
      break;
    } else {
      left_edge = i;
      sum += h_E[i];
      m++;
      mean_left = sum/m;
    }
    left = mean_left; //h_E[i];
  }

  // start from the center bin, goes to the right, look for the first drop
  // compare h_E[i] with the running mean for better edge detection

  sum = 0.0;
  m = 0;
  mean_right = right;
  for (i = center_bin+1; i < nbins; i++) {
    if (h_E[i] < (1.0-threshold)*right) {
      right_edge = i-1;
      break;
    } else {
      right_edge = i;
      sum += h_E[i];
      m++;
      mean_right = sum/m;
    }
    right = mean_right; //h_E[i];
  }

  // to be more conservative by excluding nbins_edge on both sides

  left_edge += nbins_left_edge;
  right_edge -= nbins_right_edge;
  left = h_E[left_edge];
  right = h_E[right_edge];

  mean = 0.0;
  stdev = 0.0;

  // Nbins = length(find(H));
  // NEvents = sum(H);
  int NBins = 0; 
  double NEvents = 0;
  for (i = 0; i < nbins; i++) {
    if (h_E[i] > 0) NBins++;
    NEvents += h_E[i];
  }

  d_KL = 0;
  for (i = 0; i < nbins; i++) {
    if (h_E[i] > 0) {
      double Pi = h_E[i] / NEvents;
      double Qi = 1.0 / (double)NBins;
      d_KL = d_KL + Pi*log(Pi/Qi);

      mean += h_E[i];
      stdev += h_E[i]*h_E[i];
      m++;
    }
  }

  mean /= (double)m;
  stdev = sqrt(stdev / (double)m - mean * mean);

  // is flat?
  if (d_KL < threshold)
    flat = 1;

  // print out the current histogram h_E and temperature profile T_E

  if (comm->me == 0) {
    double first_fluctuation = -1;
    write_histograms(mean, first_fluctuation, flat);
    dos_estimate();
  }

  // reduce the scaling factor
  //   if smaller than f_final, then estimate log(g(E)) from T_E
  //   otherwise, reset the energy histogram and fraction

  if (flat == 1) {
    if (f_current > f_final) {
      f_current = sqrt(f_current);
      deltaf = force->boltz * log(f_current) / (2.0 * binsize_E);
      for (i = 0; i < nbins; i++) 
        h_E[i] = 0.0;
      fraction = fraction0;
    }
  }
}

/* ---------------------------------------------------------------------- 
  Write out h_E and T_E to a text file after checking for flatness
------------------------------------------------------------------------ */

void FixSTMD::write_histograms(double mean, double fluctuation, int flat)
{
  // if one file per timestep, replace '*' with current timestep

  char *filecurrent = filename;
  char *filestar = filecurrent;
  filecurrent = new char[strlen(filestar) + 16];
  char *ptr = strchr(filestar,'*');
  *ptr = '\0';
  if (padflag == 0)
    sprintf(filecurrent,"%s" BIGINT_FORMAT "%s",
            filestar,update->ntimestep,ptr+1);
  else {
    char bif[8],pad[16];
    strcpy(bif,BIGINT_FORMAT);
    sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
    sprintf(filecurrent,pad,filestar,update->ntimestep,ptr+1);
  }
  *ptr = '*';
  
  int converged = (f_current <= f_final) ? 1 : 0;
  FILE* fp = fopen(filecurrent, "w");
  fprintf(fp, "E h_E T_E at t = " BIGINT_FORMAT " f-1 = %0.15f; "
         "edges: %d %d (%f %f); nbins = %d; mean = %f; "
         "fluctuation = %f; flat = %d; low-end flattening = %d; converged = %d; d_KL = %f\n",
         update->ntimestep, f_current-1.0, left_edge, right_edge,
         E_min + left_edge*binsize_E, E_min + right_edge*binsize_E,
         nbins, mean, fluctuation, flat, flattening, converged, d_KL);

  for (int i = 0; i < nbins; i++)
    fprintf(fp, "%g %g %g\n", E_min + i * binsize_E, h_E[i], T_E[i]);
    
  fclose(fp);
}

/* ---------------------------------------------------------------------- 
  Integrate the converged T(E) for log(g(E)) (i.e. S(E)/k_B) 
  using a finer grid of E 
------------------------------------------------------------------------ */

void FixSTMD::dos_estimate()
{
  // sanity check
  if (left_edge >= right_edge) return;

  FILE* fp = fopen("S_E.txt", "w");
  fprintf(fp, "timestep = %d; nbins = %d; scaling-1.0 = %0.15f; "
    " edges = %d %d\n", (int)update->ntimestep, 
    nbins, f_current-1.0, left_edge, right_edge);

  int i, j, iE;
  int nflatbins = right_edge-left_edge+1;
  int nebins = smooth*nflatbins;
  double dE = binsize_E / smooth;
  double left_e = E_min + left_edge * binsize_E;

  double *L = new double [nebins];
  double *Lcontin = new double [nebins];
  double *log_g = new double [nebins];
  for (i = 0; i < nebins; i++) 
    L[i] = Lcontin[i] = log_g[i] = 0.0;

  double *T = new double [nflatbins];
  for (i = 0; i < nflatbins; i++)
    T[i] = T_E[left_edge+i];

  for (iE = 0; iE < nebins; iE++) {
    j = int(iE * dE / binsize_E);
    if (j >= nflatbins-1) break;
    if (j >= 1) {
      if (T[j] != T[j-1])                  
        L[j] = binsize_E * log(T[j] / T[j-1]) / (T[j] - T[j-1]);
      else
        L[j] = binsize_E / T[j];
    }

    if (T[j+1] != T[j])
      Lcontin[iE] = (binsize_E / (T[j+1] - T[j])) * 
        log(1 + (T[j+1] - T[j]) * (iE * dE - j * binsize_E) / (binsize_E * T[j]));
    else 
      Lcontin[iE] = (iE * dE - j * binsize_E) / T[j+1];

    int idx = binsize_E * (j - 1) / dE;
    log_g[iE] = log_g[idx] + L[j] + Lcontin[iE];
    fprintf(fp, "%g %0.10f\n", iE * dE + left_e, log_g[iE]);
  }

  fclose(fp);

  double dT = 0.01;
  double T_i = T_l;

  fp = fopen("C_v.txt", "w");
  fprintf(fp, "T Cv\n");

  while (1) {

    double beta = 1/T_i;
    double F, Fmax;
    for (iE = 0; iE < nebins; iE++) {
      double e = iE * dE + left_e;
      // F = S(E) - beta E
      F = log_g[iE] - beta * e;
      if (iE == 0) {
        Fmax = F;
      } else {
        if (Fmax < F) Fmax = F;
      }
    }

    double Q = 0;
    for (iE = 0; iE < nebins; iE++) {
      double e = iE * dE + left_e;
      F = log_g[iE] - beta * e;
      Q += exp(F - Fmax);
    }

    double ave_e = 0;
    double ave_e2 = 0;
    for (iE = 0; iE < nebins; iE++) {
      double e = iE * dE + left_e;
      F = log_g[iE] - beta * e;
      ave_e += e * exp(F - Fmax) / Q;
      ave_e2 += (e*e) * exp(F - Fmax) / Q;
    }

    double Cv = (ave_e2 - ave_e*ave_e)*(beta*beta);

    fprintf(fp, "%f %f\n", T_i, Cv);

    T_i += dT;
    if (T_i > T_h) break;    
  }

  fclose(fp);

  delete [] L;
  delete [] Lcontin;
  delete [] log_g;
  delete [] T;


}

/* ----------------------------------------------------------------------
   compute energy
------------------------------------------------------------------------- */

double FixSTMD::compute_energy()
{ 
  double energy = pe->compute_scalar();
  if (ke_included) {
    double ke = temperature->compute_scalar();
    ke *= 0.5 * temperature->dof * force->boltz;
    energy += ke;
  }
  pe->addstep(update->ntimestep+1);
  if (ke_included) temperature->addstep(update->ntimestep+1);
  return energy;
} 

/* ----------------------------------------------------------------------
   resize the energy window to accomodate new values of E 
   update T_E and h_E accordingly
------------------------------------------------------------------------- */

void FixSTMD::resize_energy_window(double E)
{
   // update left or right bounds accordingly

   int j = static_cast<int>((E - E_min) / binsize_E);
   if (j >= nbins) {
     E_max += (j - nbins + 2) * binsize_E;
   } else if (j < 0) {
     E_min -= (j + 2) * binsize_E;
   }

   // grow h_E and T_E with the new number of bins

   int _nbins = ceil((E_max - E_min) / binsize_E);
   char str[256];
   sprintf(str, "Error resizing the energy window: nbins = %d; "
                "new nbins = %d; E = %f; E_min = %f; E_max = %f\n", 
           nbins, _nbins, E, E_min, E_max);
   if (_nbins <= nbins) error->all(FLERR,str);

   memory->grow(h_E, _nbins, "fix:h_e");
   memory->grow(T_E, _nbins, "fix:t_e");

   // update the resized histograms 

   if (j >= nbins) {

     // added bins are to the right

     for (int i = nbins; i < _nbins; i++) {
       h_E[i] = 0.0;
       T_E[i] = T_E[nbins-1];
     }

     if (j < _nbins && flattening == 0) 
       h_E[j] = 1.0;
   } else if (j < 0) {

     // added bins are to the left

     int naddedbins = _nbins-nbins;
     for (int i = _nbins-1; i >= 0; i--) {
       if (i >= _nbins-nbins) {
         h_E[i] = h_E[i-naddedbins];
         T_E[i] = T_E[i-naddedbins];
       } else {
         h_E[i] = 0.0;
         T_E[i] = T_E[0];
       }
     }
     int idx = (E - E_min) / binsize_E;
     if (idx >= 0 && idx < _nbins && flattening == 0) 
       h_E[idx] = 1.0;
   } 

   nbins = _nbins; 
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixSTMD::write_restart(FILE *fp)
{
  int n = 0;
  double *list = new double [2*nbins+8];
  list[n++] = T_l;
  list[n++] = T_h;
  list[n++] = E_min;
  list[n++] = E_max;
  list[n++] = binsize_E;
  list[n++] = nbins;
  list[n++] = f_current;
  list[n++] = flattening;
  for (int i = 0; i < nbins; i++) {
    list[n++] = T_E[i];
    list[n++] = h_E[i];
  }

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixSTMD::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;
  double _T_l = list[n++];
  double _T_h = list[n++];
  double _E_min = list[n++];
  double _E_max = list[n++];
  double _binsize_E = list[n++];
  int _nbins = static_cast<int> (list[n++]);
  f_current = list[n++];
  flattening = list[n++];
  if (_T_l != T_l || _T_h != T_h || _E_min != E_min || 
      _E_max != E_max || _binsize_E != binsize_E || _nbins != nbins) {
    char str[256];
    sprintf(str,"fix stmd uses parameters different from restart info: "
        "T_l = %f; T_h = %f; E_min = %f; E_max = %f; binsize = %f\n",
        _T_l, _T_h, _E_min, _E_max, _binsize_E);
    error->warning(FLERR,str);
  }

  // merge the restart data into the current
  // loop through the restart histograms
  
  int i, i_min, ibin;
  double E, T, h;
  for (i = 0; i < _nbins; i++) {

    // find the energy from the restart data

    E = _E_min + i * _binsize_E;

    // put into the current histogram, if in valid range

    ibin = (E - E_min) / binsize_E;
    T = list[n++];
    h = list[n++];

    if (ibin >= 0 && ibin < nbins) {
      T_E[ibin] = T;
      h_E[ibin] = h;
      if (i == 0) i_min = ibin;
    }
  }

  // flatten the low-energy end
  
  for (i = 0; i < i_min; i++) {
    T_E[i] = T_E[i_min];
  }
}

/* ---------------------------------------------------------------------- */

int FixSTMD::modify_param(int narg, char **arg)
{
  int iarg=0;
  if (strcmp(arg[iarg],"fd") == 0) {
    if (narg < 3) error->all(FLERR,"Illegal fix_modify command");
    double fd = force->numeric(FLERR,arg[iarg+1]);
    if (fd < 0.0) error->all(FLERR,"fix stmd: initial fd should be positive");
    f_current = 1.0 + fd;
    fd = force->numeric(FLERR,arg[iarg+2]);
    if (fd < 0.0) error->all(FLERR,"fix stmd: final fd should be positive");
    f_final = 1.0 + fd;
    return 3;

  } else if (strcmp(arg[iarg],"flatcheck") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    flatcheck = force->numeric(FLERR,arg[iarg+1]);
    return 2;

  } else if (strcmp(arg[iarg],"flattening") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    flattening = force->numeric(FLERR,arg[iarg+1]);
    return 2;

  } else if (strcmp(arg[iarg],"left_edge") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    nbins_left_edge = force->numeric(FLERR,arg[iarg+1]);
    return 2;

  } else if (strcmp(arg[iarg],"right_edge") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    nbins_right_edge = force->numeric(FLERR,arg[iarg+1]);
    return 2;

  } else if (strcmp(arg[iarg],"threshold") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    threshold = force->numeric(FLERR,arg[iarg+1]);
    return 2;

  } else if (strcmp(arg[iarg],"smooth") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    smooth = force->numeric(FLERR,arg[iarg+1]);
    return 2;

  } else if (strcmp(arg[iarg],"fraction") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    fraction0 = force->numeric(FLERR,arg[iarg+1]);
    fraction = fraction0;
    return 2;

  } else if (strcmp(arg[iarg],"maxcount") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    maxcount = force->numeric(FLERR,arg[iarg+1]);
    return 2;
 
  } else if (strcmp(arg[iarg],"file") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    delete [] filename;
    int n = strlen(arg[iarg+1]) + 1;
    filename = new char[n];
    strcpy(filename,arg[iarg+1]);
    return 2;

  } else if (strcmp(arg[iarg],"pad") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    padflag = force->inumeric(FLERR,arg[iarg+1]);
    if (padflag < 0) error->all(FLERR,"Illegal fix stmd command");
    return 2;

  } else if (strcmp(arg[iarg],"fixed") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    fixed = force->inumeric(FLERR,arg[iarg+1]);
    return 2;
  } 
  return 0;
}

/* ----------------------------------------------------------------------
   current scaling factor less 1, i.e., f_d
------------------------------------------------------------------------- */

double FixSTMD::compute_scalar()
{
  return (f_current-1.0);
}

/* ----------------------------------------------------------------------
   return corresponding parameters
     n = 0: T_current
     n = 1: E_current
     n = 2: 1/0 for convergence or not 
     n = 3: current scaling factor
     n = 4: nout
   to print out these parameters by thermo, use (n+1) as index
   For example, f_stmd[2] corresponds to n = 1
------------------------------------------------------------------------- */

double FixSTMD::compute_vector(int n)
{
  double value;
  if (n == 0) value = T_current;
  else if (n == 1) value = E_current;
  else if (n == 2) value = (double)((f_current <= f_final) ? 1 : 0);
  else if (n == 3) value = f_current;
  else if (n == 4) value = nout;
  else value = -1.0;
  return value;
}

/* ----------------------------------------------------------------------
   The below implementation is closer to what's described in 
   Refs. Kim et al. Phys. Rev. Lett. 97, 050601 (2006)
     and Kim et al. J. Chem. Phys. 126, 135101 (2007)
   but does not take into account limiting cases (e.g. T[j+1] == T[j])
------------------------------------------------------------------------- */
/* 
void FixSTMD::dos_estimate()
{
  // sanity check
  if (left_edge >= right_edge) return;

  FILE* fp = fopen("S_E.txt", "w");
  fprintf(fp, "timestep = %d; nbins = %d; scaling-1.0 = %lf; "
    " edges = %d %d\n", (int)update->ntimestep, 
    nbins, f_current-1.0, left_edge, right_edge);

  int i, j, ibin, istar, nebins;
  double E, Ei, Eim1, Eip1, Ebim1, Ebi, dE;

  // find low energy bound at the left edge 

  double E_l = E_min + left_edge * binsize_E;

  // use a finer bin size

  nebins = smooth*nbins;
  dE = binsize_E / smooth;

  double *S = new double [nebins];   
  int l = left_edge;

  // L and lambda are evaluted at each grid points of 
  // the original energy histogram
  double *L = new double [nbins];
  double *lambda = new double [nbins];
  for (j = 0; j < nbins; j++) 
    lambda[j] = L[j] = 0.0;
  for (j = left_edge; j <= right_edge; j++) {
    lambda[j] = (T_E[j+1] - T_E[j]) / binsize_E;
    if (j > left_edge)
      L[j] = log(1 + lambda[j-1]*binsize_E/T_E[j-1]) / lambda[j-1];
  }

  for (int n = 0; n<nebins; n++) {

    // energy at the finer grid for interpolation

    E = E_l + n * dE;

    // energy at the coarse grid used for STMD

    i = (E - E_min) / binsize_E;

    // E[i-1], E[i] and E[i+1]

    Eim1 = E_min + (i - 1) * binsize_E;
    Ei   = E_min + i * binsize_E;
    Eip1 = E_min + (i + 1) * binsize_E;

    // middle values

    Ebim1 = 0.5 * (Eim1 + Ei);
    Ebi = 0.5 * (Ei + Eip1);

    // find the value of i^*
    if (E >= Ebim1 && E <= Ei) istar = i -1;
    if (E >= Ei && E <= Ebi) istar = i;

    // accumulate to S[n]
    S[n] = log(1 + lambda[i]*(E - Ei)/T_E[i]) / lambda[i];
    for (j = l+1; j <= istar; j++)       
      S[n] += L[j];
      
    fprintf(fp, "%f %f\n", E, S[n]);
  }

  delete [] S;
  delete [] lambda;
  delete [] L;

  fclose(fp);
}
*/


