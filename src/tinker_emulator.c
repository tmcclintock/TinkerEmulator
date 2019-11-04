#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//Squared exponential kernel function
double kernel(double*x1, double*x2, double*lambdas){
  return 0;
}

//Emulate the sub-parameters; void for now
//double* emulate(double*cos, double*subpars){
void emulate(double*cosmological_parameters, double*subpars){
  int i, j;
  return;
}

//Predict tinker parameters; void for now
//double* predict_tinker_parameters(double*cos, double redshift, double*tpars){
//void predict_tinker_parameters(double*cosmological_parameters,
//			       double redshift, double*tinker_parameters){
//  int i,j;
//  return;
//}

typedef struct{
    //parameters of equation 7 in 1907.13167.
    double A; 
    double a; 
    double B; 
    double b;
    double C;
    double c;
} Emu_tinker_param; 

typedef struct{
    double ombh2;  
    double omch2; 
    double w0; 
    double n_s;
    double H0; 
    double ln_10_10_A_s;  //ln(10^10As)
    double N_eff;
} Emu_tinker_cosmo ;
//Predict tinker parameters;
Emu_tinker_param predict_tinker_parameters(Emu_tinker_cosmo emutinker, double redshift){
    Emu_tinker_param test;
    return test; 
}
int main(){
  printf("test\n");
  return 0;
}
