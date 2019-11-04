"""Emulates parameters in the Tinker mass function.
"""
import numpy as np
import scipy as sp
import george

class TinkerEmulator(object):

    def __init__(self):
        self.pars = np.loadtxt("./data/parameters.txt")
        self.means = np.loadtxt("./data/rotated_means.txt")
        self.variances = np.loadtxt("./data/rotated_variances.txt")
        self.rotation_matrix = np.loadtxt("./data/rotation_matrix.txt")

        #attributes that will be filled
        self.GP_list = []
        
        #write this function
        self._train()

    def _train(self):
        hps = np.std(self.pars, 0) #Guess for hyperparameters
        means = self.means
        stds = np.sqrt(self.variances)
        N = len(means[0])
        Ncos = len(hps)
        
        for i in range(N):
            k = george.kernels.ExpSquaredKernel(hps, ndim=Ncos)
            gp = george.GP(k, mean=np.mean(means[:,i]))
            gp.compute(self.pars, stds[:,i])
            self.GP_list.append(gp)
        return
        
    def emulate(self, cosmological_parameters):
        """Emulate the sub-parameters in the Tinker mass function as
        implemented in the McClintock+ (2019) halo bias emulator.

        TODO: add the option to return the uncertainties

        Args:
            cosmological_parameters (array-like): 
                Obh2 Och2 w0 n_s ln10As H0 N_eff

        Returns:
            B0 c0 A1 B1

        """
        means = self.means.T
        R = self.rotation_matrix
        p_r = np.array([gp.predict(y, cosmological_parameters)[0]
                           for y,gp in zip(means, self.GP_list)])
        return np.dot(R, p_r).flatten()

    def predict_tinker_parameters(self, cosmological_parameters, redshift):
        """Predict the Tinker parameters.

        Args:
            cosmological_parameters (array-like): 
                Obh2 Och2 w0 n_s ln10As H0 N_eff
            redshift (float): cosmological redshift `z`

        Returns:
            A a B b C c tinker parameters

        """
        x = 1. / (1. + redshift) - 0.5
        A0 = 4.2828605
        a0 = 0.4722138
        b0 = 1.5170196
        C0 = 0.888452
        a1 = -0.56318698
        b1 = -0.63010135
        C1 = -0.5956625
        c1 = -1.85148405
        B0, c0, A1, B1 = self.emulate(cosmological_parameters)
        A = A0 + x * A1
        a = a0 + x * a1
        B = B0 + x * B1
        b = b0 + x * b1
        C = C0 + x * C1
        c = c0 + x * c1
        return A,a,B,b,C,c
