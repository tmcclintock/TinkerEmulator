//Squared exponential kernel function
double kernel(double *x1, double *x2, double *loglambdas, int ndim){
    double k=0.0;
    int i;
    for(i=0; i<ndim; ++i){ 
        k += pow(x1[i]- x2[i], 2)/(2*exp(loglambdas[i]));
    }
    return exp(-1.0*k);
}
void create_invese_covariance_matrix(double **  sigma_x_x_inv, int nparam){
    int i,j;
    double result;
    gsl_matrix * cov   = gsl_matrix_calloc(tinkerEmuParam.tinker_bias_nsamp, tinkerEmuParam.tinker_bias_nsamp);
    gsl_matrix_set_zero(cov);
    for(i=0; i<tinkerEmuParam.tinker_bias_nsamp; ++i){
        for(j=0; j<tinkerEmuParam.tinker_bias_nsamp; ++j){
            result = kernel(emu_tinker_bias_cosmopars[i], emu_tinker_bias_cosmopars[j], emu_tinker_bias_log_lambda[nparam], tinkerEmuParam.tinker_bias_ncosmo); 
            if(i==j) result += emu_tinker_bias_variances[i][nparam];
            gsl_matrix_set(cov,i,j,result);
        } 
    }
    invert_matrix_colesky(cov);
    for(i=0; i<tinkerEmuParam.tinker_bias_nsamp; ++i){
        for(j=0; j<tinkerEmuParam.tinker_bias_nsamp; ++j){
            sigma_x_x_inv[i][j] = gsl_matrix_get(cov,i,j);
        }
    } 
    gsl_matrix_free(cov);
}

double bias_emulate_single_param(double*cos, int nparam){
    int i,j;
    static double ** sigma_x_x_inv0  = 0;
    static double ** sigma_x_x_inv1  = 0;
    static double ** sigma_x_x_inv2  = 0;
    static double ** sigma_x_x_inv3  = 0;
    if(sigma_x_x_inv0==0){
        sigma_x_x_inv0=create_double_matrix(0, tinkerEmuParam.tinker_bias_nsamp, 0, tinkerEmuParam.tinker_bias_nsamp);
        sigma_x_x_inv1=create_double_matrix(0, tinkerEmuParam.tinker_bias_nsamp, 0, tinkerEmuParam.tinker_bias_nsamp);
        sigma_x_x_inv2=create_double_matrix(0, tinkerEmuParam.tinker_bias_nsamp, 0, tinkerEmuParam.tinker_bias_nsamp);
        sigma_x_x_inv3=create_double_matrix(0, tinkerEmuParam.tinker_bias_nsamp, 0, tinkerEmuParam.tinker_bias_nsamp);
        create_invese_covariance_matrix(sigma_x_x_inv0, 0);
        create_invese_covariance_matrix(sigma_x_x_inv1, 1);
        create_invese_covariance_matrix(sigma_x_x_inv2, 2);
        create_invese_covariance_matrix(sigma_x_x_inv3, 3);
    }
    double ** sigma_x_x_inv=0;
    switch(nparam){
        case 0: sigma_x_x_inv = sigma_x_x_inv0; break;
        case 1: sigma_x_x_inv = sigma_x_x_inv1; break;
        case 2: sigma_x_x_inv = sigma_x_x_inv2; break;
        case 3: sigma_x_x_inv = sigma_x_x_inv3; break;
        default:
            printf("bias_emulate_single_param: we only have %d parameters but nparam: %d is passed", tinkerEmuParam.tinker_bias_nparam, nparam);
    }
    if(nparam==0){
    }
    double result=0;
    double temp = 0;
    for(i=0; i<tinkerEmuParam.tinker_bias_nsamp; ++i){
        temp = 0;
        for(j=0; j<tinkerEmuParam.tinker_bias_nsamp; ++j){
           temp += (sigma_x_x_inv)[i][j]*emu_tinker_bias_y[j][nparam]; 
        }
        temp *= kernel(cos, emu_tinker_bias_cosmopars[i], emu_tinker_bias_log_lambda[nparam], tinkerEmuParam.tinker_bias_ncosmo);
        result += temp ;
    }
    result += emu_tinker_bias_ymean[nparam];
    return result;
}

void bias_emulate_all(double*cos, double*ystar){
    int i,j;
    double * yin = create_double_vector(0, tinkerEmuParam.tinker_bias_nparam);
    for(i=0; i<tinkerEmuParam.tinker_bias_nparam; ++i){
        yin[i] = bias_emulate_single_param(cos, i);
    }
    //rotate:
    for(i=0; i<tinkerEmuParam.tinker_bias_nparam; ++i){
        for(j=0; j<tinkerEmuParam.tinker_bias_nparam; ++j){
            ystar[i] += emu_tinker_bias_rotation_matrix[i][j]*yin[j];
        }
    }   
      
    
    
    free_double_vector(yin, 0, tinkerEmuParam.tinker_bias_nparam);
}


typedef struct{
    double A;
    double a;
    double B;
    double b;
    double C;
    double c;
    //parameters of equation 7 in 1907.13167.;
} Emu_tinker_bias_param; 

//Predict tinker parameters;
void predict_tinker_bias_parameters(double ombh2, double omch2, double w0, double n_s, double H0, double ln_10_10_A_s, double N_eff,  double redshift, double* emu_tinker_bias_param){
    double cosmo_emu[7] = {ombh2, omch2, w0, n_s, H0, ln_10_10_A_s, N_eff};
    double *ystar = create_double_vector(0, tinkerEmuParam.tinker_bias_nparam);
    bias_emulate_all(cosmo_emu, ystar);
    double x = 1. / (1. + redshift) - 0.5;
    double A0 = 4.2828605;
    double a0 = 0.4722138;
    double b0 = 1.5170196;
    double C0 = 0.888452;
    double a1 = -0.56318698;
    double b1 = -0.63010135;
    double C1 = -0.5956625;
    double c1 = -1.85148405;
    emu_tinker_bias_param[0] = A0 + x * ystar[2];
    emu_tinker_bias_param[1] = a0 + x * a1;
    emu_tinker_bias_param[2]= ystar[0] + x * ystar[3];
    emu_tinker_bias_param[3] = b0 + x * b1;
    emu_tinker_bias_param[4]= C0 + x * C1;
    emu_tinker_bias_param[5] = ystar[1] + x * c1;
    free_double_vector(ystar, 0, tinkerEmuParam.tinker_bias_nparam);
}
