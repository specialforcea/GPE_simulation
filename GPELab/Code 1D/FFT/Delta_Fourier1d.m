%% Computation of the Fourier frequencies for the Laplacian
%% INPUT:
%%          Geometry1D: Structure containing variables concerning the geometry of the problem in 1D in the FFT context (structure) (see FFTGeometry1D_Var1d.m)
%% OUTPUT:
%%          Dx frequencies for the Laplacian in the x direction

function [Dx]=Delta_Fourier1d(FFTGeometry1D)
%% Computation of the frequencies for the gradient
[Gx] = (2*pi*1i.*[[0:ceil(FFTGeometry1D.Nx/2)-1],[-floor(FFTGeometry1D.Nx/2):1:-1]]/FFTGeometry1D.Lx)';
Dx = Gx.^2;