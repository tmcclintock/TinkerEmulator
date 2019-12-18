//Squared exponential kernel function
double kernel(double *x1, double *x2, double *loglambdas, int ndim){
    double k=0.0;
    int i;
    for(i=0; i<ndim; ++i){ 
        k += pow(x1[i]- x2[i], 2)/(2*exp(loglambdas[i]));
    }
    return exp(-1.0*k);
}
void create_invese_covariance_matrix(double **  sigma_x_x_inv, int nparam, int type){
    int i,j, k;
    double result;
    double (*cosmopars)[7];
    double (*loglambda)[7];
    int nsamp, ncosmo;
    double mean_std = 0;
    switch(type){
        case 0: //bias emulator 
            cosmopars = emu_tinker_bias_cosmopars; 
            loglambda = emu_tinker_bias_log_lambda;
            nsamp = tinkerEmuParam.tinker_bias_nsamp; 
            ncosmo = tinkerEmuParam.tinker_bias_ncosmo;
            break;
        case 1://hmf emulator
            cosmopars = emu_tinker_hmf_cosmopars; 
            loglambda = emu_tinker_hmf_log_lambda;
            nsamp = tinkerEmuParam.tinker_hmf_nsamp; 
            ncosmo = tinkerEmuParam.tinker_hmf_ncosmo;
            break;
        default: 
            printf("type value %d is passed to create_invese_covariance_matrix, but we only support 0 [bias emulator] or 1 [hmf emulator]\n");
            exit(0);
    }
    gsl_matrix * cov   = gsl_matrix_calloc(nsamp, nsamp);
    gsl_matrix_set_zero(cov);
    for(i=0; i<nsamp; ++i){
        for(j=0; j<nsamp; ++j){
            result = kernel(cosmopars[i], cosmopars[j], loglambda[nparam], ncosmo); 
            if(i==j){
                if(type==0) result += emu_tinker_bias_variances[i][nparam];
                if(type==1){
                //this is for white_noise agrument in the hmf emulator
                    result += emu_tinker_hmf_variances[i][nparam];
                    mean_std = 0; 
                    for(k=0; k<nsamp; ++k){
                        mean_std += sqrt(emu_tinker_hmf_variances[k][nparam]); 
                    } 
                    mean_std /= (double) nsamp;
                    result += pow(mean_std, 2);
                }
            }
            gsl_matrix_set(cov,i,j,result);
        } 
    }
    invert_matrix_colesky(cov);
    for(i=0; i<nsamp; ++i){
        for(j=0; j<nsamp; ++j){
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
        create_invese_covariance_matrix(sigma_x_x_inv0, 0, 0);
        create_invese_covariance_matrix(sigma_x_x_inv1, 1, 0);
        create_invese_covariance_matrix(sigma_x_x_inv2, 2, 0);
        create_invese_covariance_matrix(sigma_x_x_inv3, 3, 0);
    }
    double ** sigma_x_x_inv=0;
    switch(nparam){
        case 0: sigma_x_x_inv = sigma_x_x_inv0; break;
        case 1: sigma_x_x_inv = sigma_x_x_inv1; break;
        case 2: sigma_x_x_inv = sigma_x_x_inv2; break;
        case 3: sigma_x_x_inv = sigma_x_x_inv3; break;
        default:
            printf("bias_emulate_single_param: we only have %d parameters but nparam: %d is passed\n", tinkerEmuParam.tinker_bias_nparam, nparam);
            exit(0);
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


double hmf_emulate_single_param(double*cos, int nparam){
    int i,j;
    static double ** sigma_x_x_inv0  = 0;
    static double ** sigma_x_x_inv1  = 0;
    static double ** sigma_x_x_inv2  = 0;
    static double ** sigma_x_x_inv3  = 0;
    static double ** sigma_x_x_inv4  = 0;
    static double ** sigma_x_x_inv5  = 0;
    if(sigma_x_x_inv0==0){
        sigma_x_x_inv0=create_double_matrix(0, tinkerEmuParam.tinker_hmf_nsamp, 0, tinkerEmuParam.tinker_hmf_nsamp);
        sigma_x_x_inv1=create_double_matrix(0, tinkerEmuParam.tinker_hmf_nsamp, 0, tinkerEmuParam.tinker_hmf_nsamp);
        sigma_x_x_inv2=create_double_matrix(0, tinkerEmuParam.tinker_hmf_nsamp, 0, tinkerEmuParam.tinker_hmf_nsamp);
        sigma_x_x_inv3=create_double_matrix(0, tinkerEmuParam.tinker_hmf_nsamp, 0, tinkerEmuParam.tinker_hmf_nsamp);
        sigma_x_x_inv4=create_double_matrix(0, tinkerEmuParam.tinker_hmf_nsamp, 0, tinkerEmuParam.tinker_hmf_nsamp);
        sigma_x_x_inv5=create_double_matrix(0, tinkerEmuParam.tinker_hmf_nsamp, 0, tinkerEmuParam.tinker_hmf_nsamp);
        create_invese_covariance_matrix(sigma_x_x_inv0, 0, 1);
        create_invese_covariance_matrix(sigma_x_x_inv1, 1, 1);
        create_invese_covariance_matrix(sigma_x_x_inv2, 2, 1);
        create_invese_covariance_matrix(sigma_x_x_inv3, 3, 1);
        create_invese_covariance_matrix(sigma_x_x_inv4, 4, 1);
        create_invese_covariance_matrix(sigma_x_x_inv5, 5, 1);
    }
    double ** sigma_x_x_inv=0;
    switch(nparam){
        case 0: sigma_x_x_inv = sigma_x_x_inv0; break;
        case 1: sigma_x_x_inv = sigma_x_x_inv1; break;
        case 2: sigma_x_x_inv = sigma_x_x_inv2; break;
        case 3: sigma_x_x_inv = sigma_x_x_inv3; break;
        case 4: sigma_x_x_inv = sigma_x_x_inv4; break;
        case 5: sigma_x_x_inv = sigma_x_x_inv5; break;
        default:
            printf("hmf_emulate_single_param: we only have %d parameters but nparam: %d is passed\n", tinkerEmuParam.tinker_hmf_nparam, nparam);
            exit(0);
    }
    double result=0;
    double temp = 0;
    for(i=0; i<tinkerEmuParam.tinker_hmf_nsamp; ++i){
        temp = 0;
        for(j=0; j<tinkerEmuParam.tinker_hmf_nsamp; ++j){
           temp += (sigma_x_x_inv)[i][j]*emu_tinker_hmf_y[j][nparam]; 
        }
        temp *= kernel(cos, emu_tinker_hmf_cosmopars[i], emu_tinker_hmf_log_lambda[nparam], tinkerEmuParam.tinker_hmf_ncosmo);
        result += temp ;
    }
    result += emu_tinker_hmf_ymean[nparam];
    return result;
}


void emulate_all(double*cos, double*ystar, int type){
    int i,j;
    int nparam;
    switch(type){
        case 0: //bias emulator 
            nparam = tinkerEmuParam.tinker_bias_nparam;
            break;
        case 1://hmf emulator
            nparam = tinkerEmuParam.tinker_hmf_nparam;
            break;
        default: 
            printf("type value %d is passed to create_invese_covariance_matrix, but we only support 0 [bias emulator] or 1 [hmf emulator]\n");
            exit(0);
    }
    
    double * yin = create_double_vector(0, nparam);
    for(i=0; i<nparam; ++i){
        if(type==0) yin[i] = bias_emulate_single_param(cos, i);
        if(type==1) yin[i] = hmf_emulate_single_param(cos, i);
    }
    //rotate:
    for(i=0; i<nparam; ++i){
        for(j=0; j<nparam; ++j){
            if(type==0) ystar[i] += emu_tinker_bias_rotation_matrix[i][j]*yin[j];
            if(type==1) ystar[i] += emu_tinker_hmf_rotation_matrix[i][j]*yin[j];
        }
    }   
    free_double_vector(yin, 0, nparam);
}


//Predict tinker parameters;
void predict_tinker_bias_parameters(double ombh2, double omch2, double w0, double n_s, double H0, double ln_10_10_A_s, double N_eff,  double redshift, double* emu_tinker_bias_param){
    double cosmo_emu[7] = {ombh2, omch2, w0, n_s, H0, ln_10_10_A_s, N_eff};
    double *ystar = create_double_vector(0, tinkerEmuParam.tinker_bias_nparam);
    emulate_all(cosmo_emu, ystar, 0);
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

//Predict tinker hmf parameters;
void predict_tinker_hmf_parameters(double ombh2, double omch2, double w0, double n_s, double H0, double ln_10_10_A_s, double N_eff,  double redshift, double* emu_tinker_hmf_param){
    double cosmo_emu[7] = {ombh2, omch2, w0, n_s, H0, ln_10_10_A_s, N_eff};
    double *ystar = create_double_vector(0, tinkerEmuParam.tinker_hmf_nparam);//e0,f0,g0,d1,e1,g1
    emulate_all(cosmo_emu, ystar, 1);
    double x = 1. / (1. + redshift) - 0.5;
    double d0 = 2.39279115;
    double f1 =  0.11628991;
    emu_tinker_hmf_param[0] = d0 + x * ystar[3];        //d
    emu_tinker_hmf_param[1] = ystar[0] + x * ystar[4];  //e
    emu_tinker_hmf_param[2] = ystar[1] + x * f1;        //f
    emu_tinker_hmf_param[3] = ystar[2] + x * ystar[5];  //g
    free_double_vector(ystar, 0, tinkerEmuParam.tinker_hmf_nparam);
}






