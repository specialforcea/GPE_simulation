%% Computation of the Fourier frequencies for the gradient
%% INPUT:
%%          Geometry2D: Structure containing variables concerning the geometry of the problem in 2D in the FFT context (structure) (see FFTGeometry2D_Var2d.m)
%% OUTPUT:
%%          Gx,Gy: frequencies for the gradient in the x,y direction

function [Gx,Gy]=Grad_Fourier2d(FFTGeometry2D)
%% Computation of the frequencies for the gradient
[Gx,Gy] = meshgrid((2*pi*1i.*[[0:ceil(FFTGeometry2D.Nx/2)-1],[-floor(FFTGeometry2D.Nx/2):1:-1]]/FFTGeometry2D.Lx),(2*pi*1i.*[[0:ceil(FFTGeometry2D.Ny/2)-1],[-floor(FFTGeometry2D.Ny/2):1:-1]]/FFTGeometry2D.Ly));