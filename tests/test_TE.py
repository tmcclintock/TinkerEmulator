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

