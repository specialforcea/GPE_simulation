%% Computation of the Fourier frequencies for the Laplacian
%% INPUT:
%%          Geometry2D: Structure containing variables concerning the geometry of the problem in 2D in the FFT context (structure) (see FFTGeometry2D_Var2d.m)
%% OUTPUT:
%%          Dx,Dy: frequencies for the Laplacian in the x,y direction

function [Dx,Dy]=Delta_Fourier2d(FFTGeometry2D)
%% Computation of the frequencies for the gradient
[Gx,Gy] = meshgrid((2*pi*1i.*[[0:ceil(FFTGeometry2D.Nx/2)-1],[-floor(FFTGeometry2D.Nx/2):1:-1]]/FFTGeometry2D.Lx),(2*pi*1i.*[[0:ceil(FFTGeometry2D.Ny/2)-1],[-floor(FFTGeometry2D.Ny/2):1:-1]]/FFTGeometry2D.Ly));
Dx = Gx.^2;
Dy = Gy.^2;