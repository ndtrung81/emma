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

#ifdef FIX_CLASS

FixStyle(stmd,FixSTMD)

#else

#ifndef LMP_FIX_STMD_H
#define LMP_FIX_STMD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSTMD : public Fix {
 public:
  FixSTMD(class LAMMPS *, int, char **);
  ~FixSTMD();
  int setmask();
  void init();
  void setup(int);
  void post_force(int);
  void write_restart(FILE *fp);
  void restart(char *buf);
  int modify_param(int, char **);
  double compute_scalar();
  double compute_vector(int n);

 private:
  int nbins;                             // number of energy bins
  double T_l, T_h;                       // temperature window
  double T_0;                            // set temperature used for maintaining the kinetic energy
  double T_current;                      // configuration temperature interpolated from T_E at current potential energy
  double *T_E;                           // temperature estimate (i.e. temperature histogram)
  double binsize_E;                      // energy bin size
  double E_min, E_max;                   // energy window
  double E_current;                      // current energy
  double *h_E;                           // energy histogram
  double f_current, f_final, deltaf;     // scaling factors and deltaf = log(f)/(2*binsizeE)
  double threshold;                      // fraction of the mean value, below which the histogram is considered flat
  double fraction, fraction0;            // fraction of nbins between left and right edges to check flatness
  int maxcount;                          // maximum value of the average of the energy bins to adjust fraction
  int flatcheck;                         // period to check flatness
  int flattening;                        // 1/0 to flattening the low-energy end of T_E (during the inital stage)
  int nbins_left_edge, nbins_right_edge; // number of excluded bins at two ends after flatness check
  int left_edge, right_edge;             // two edges to check flatness and integrate for the DOS
  int smooth;                            // finer resolution for the energy grid used for integration for the DOS
  char* filename;                        // file name to write out histograms
  int padflag;                           // number of characters to convert the time step in the file name to
  int fixed;                             // 1 if the energy window is fixed, 0 otherwise
  int ke_included;                       // 1 if kinetic energy is included to the energy, 0 otherwise
  int first_time;                        // 1 if the energy is in the window for the first time, 0 otherwise
  int nout;                              // number of steps out of the energy window
  double d_KL;                           // Kullbach-Leibler flatness metrics

  char *id_pe;                           
  class Compute *pe;
  char *id_temp;
  class Compute *temperature;
  
  int check_flatness();
  int check_flatness_kl();
  void dos_estimate();
  double compute_energy();
  void low_end_flattening();
  void resize_energy_window(double E);
  void write_histograms(double mean, double fluctuation, int flat);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot use fix enforce2d with 3d simulation

Self-explanatory.

*/
