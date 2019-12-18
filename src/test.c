#include <math.h>
#include <stdlib.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <fftw3.h>


#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_deriv.h>

#include "../../../cosmolike_core/theory/basics.c"
#include "tinker_emulator.h"
#include "tinker_emulator.c"

int main(){
  //Test kernel
  double x0[7] = {0.0043*pow(0.7,2), 0.286*pow(0.7,2), -1, 1, 3.0, 70, 3.0};
  double x1[7] = {0.0043*pow(0.7,2), 0.286*pow(0.7,2), -1, 1, 3.0, 70, 5.0};
  printf("kernel : %e \n", kernel(x0, x1, emu_tinker_bias_log_lambda[0], tinkerEmuParam.tinker_bias_ncosmo));
  //Test bias emulator 
  double emu_tinker_bias_param[6];
  predict_tinker_bias_parameters(0.0043*pow(0.7,2), 0.286*pow(0.7,2), -1, 1, 3.0, 70, 3.0, 0.3, emu_tinker_bias_param);
  double bias_ref[6] = {4.432419358797392, 0.32058653615384625, 1.3586846172228393, 1.347376928846154, 0.728081326923077, -1.0947315110189246}; 
  printf("bias result:\n");
  printf("item c ref\n");
  for(int i=0; i<tinkerEmuParam.tinker_bias_nparam_redshift; ++i) printf("%d, %e, %e\n", i,emu_tinker_bias_param[i], bias_ref[i]);
  printf("####################\n");
  
  //Test hmf emulator
  double emu_tinker_hmf_param[4];
  double hmf_ref[4] = { 2.5198116919192444, 0.8454997909106351, 0.5441108055593771, 1.3001900391246013}; 
  predict_tinker_hmf_parameters(0.0043*pow(0.7,2), 0.286*pow(0.7,2), -1, 1, 3.0, 70, 3.0, 0.3, emu_tinker_hmf_param);
  printf("hmf result:\n");
  printf("item c ref\n");
  for(int i=0; i<tinkerEmuParam.tinker_hmf_nparam_redshift; ++i) printf("%d, %e, %e\n", i,emu_tinker_hmf_param[i], hmf_ref[i]);
  printf("####################\n");
  
    
}
