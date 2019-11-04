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

Emu_tinker_param predict_tinker_parameters(Emu_tinker_cosmo, double redshift){

}
