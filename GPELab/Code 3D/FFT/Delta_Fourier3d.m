%% Computation of the Fourier frequencies for the Laplacian
%% INPUT:
%%          Geometry3D: Structure containing variables concerning the geometry of the problem in 3D in the FFT context (structure) (see FFTGeometry3D_Var3d.m)
%% OUTPUT:
%%          Dx,Dy,Dz: frequencies for the Laplacian in the x,y,z direction

function [Dx,Dy,Dz]=Delta_Fourier3d(FFTGeometry3D)
%% Computation of the frequencies for the gradient
[Gx,Gy,Gz] = meshgrid((2*pi*1i.*[[0:ceil(FFTGeometry3D.Nx/2)-1],[-floor(FFTGeometry3D.Nx/2):1:-1]]/FFTGeometry3D.Lx),(2*pi*1i.*[[0:ceil(FFTGeometry3D.Ny/2)-1],[-floor(FFTGeometry3D.Ny/2):1:-1]]/FFTGeometry3D.Ly),(2*pi*1i.*[[0:ceil(FFTGeometry3D.Nz/2)-1],[-floor(FFTGeometry3D.Nz/2):1:-1]]/FFTGeometry3D.Lz));
Dx = Gx.^2;
Dy = Gy.^2;
Dz = Gz.^2;