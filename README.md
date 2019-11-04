# TinkerEmulator
Emulator for the Tinker halo bias parameters, without involving a Boltzmann code.

This package will support the stand-alone parameter emulator in both Python and C (since it will be incorporated into CosmoLike). Here is a list of TODO items:

1) Copy over the data and parameter files from the main emulator repository.
2) Build an emulator in Python with `george`. In a notebook, print and save the kernel parameters.
3) Implement the Gaussian process in C and built the emulator.