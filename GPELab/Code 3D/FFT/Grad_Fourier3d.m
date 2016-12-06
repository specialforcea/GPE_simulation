%% Computation of the Fourier frequencies for the gradient
%% INPUT:
%%          Geometry3D: Structure containing variables concerning the geometry of the problem in 3D in the FFT context (structure) (see FFTGeometry3D_Var3d.m)
%% OUTPUT:
%%          Gx,Gy,Gz: frequencies for the gradient in the x,y,z direction

function [Gx,Gy,Gz]=Grad_Fourier3d(FFTGeometry3D)
%% Computation of the frequencies for the gradient
[Gx,Gy,Gz] = meshgrid((2*pi*1i.*[[0:ceil(FFTGeometry3D.Nx/2)-1],[-floor(FFTGeometry3D.Nx/2):1:-1]]/FFTGeometry3D.Lx),(2*pi*1i.*[[0:ceil(FFTGeometry3D.Ny/2)-1],[-floor(FFTGeometry3D.Ny/2):1:-1]]/FFTGeometry3D.Ly),(2*pi*1i.*[[0:ceil(FFTGeometry3D.Nz/2)-1],[-floor(FFTGeometry3D.Nz/2):1:-1]]/FFTGeometry3D.Lz));