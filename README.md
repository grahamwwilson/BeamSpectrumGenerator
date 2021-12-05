# BeamSpectrumGenerator
Generate (x1,x2) distributions

## test.f
Initial approach of test.f simply uses 3-parameter CIRCE1 F77 implementation 
of the beta distribution generator for sampling (x1,x2) from a CIRCE function 
that describes the scaled beam energy after potential beamstrahlung.

## testb.f
testb.f produces (x1,x2) after convolving 
Gaussian beam-energy-spread (BES) with beamstrahlung.
Again uses 3-parameter CIRCE1 F77 implementation of the beta distribution 
generator for sampling from the CIRCE function for the beamstrahlung part.
