"""Tests of the Tinker emulator.
"""
import TinkerEmulator
import numpy as np
import scipy as sp
import numpy.testing as npt
import pytest

def test_TE_basic():
    #Smoke tests
    TE = TinkerEmulator.TinkerEmulator()

def test_TE_prediction():
    #Mean of all cosmological parameters
    p = np.array([ 2.22629225e-02,  1.17830325e-01, -9.96512725e-01,  9.62515150e-01,   3.08894925e+00,  6.82317125e+01,  3.45000000e+00])

    TE = TinkerEmulator.TinkerEmulator()
    #Test that the prediction is being made
    ps = TE.emulate(p)
    npt.assert_equal(len(ps), 4)

    zs = [0, 0.5, 1., 2]
    for z in zs:
        tps = TE.predict_tinker_parameters(p, z)
        npt.assert_equal(len(tps), 6)
    
if __name__ == "__main__":
    test_TE_basic()
    test_TE_prediction()
