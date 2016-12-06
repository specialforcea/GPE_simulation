%% Creation of the 3D geometry structure for the FFT
%% INPUT:
%%          Geometry3D: Structure containing variables concerning the 3D geometry (structure) (see Geometry3D_Var3d.m)
%% OUTPUT:
%%          FFTGeometry3D: Structure containing variables concerning the 3D geometry for the FFT (structure)

function [FFTGeometry3D] = FFTGeometry3D_Var3d(Geometry3D)
%% Initialization
FFTGeometry3D = Geometry3D; % Copying the 3D geometry

%% Preparation for the FFT of the 3D geometry
FFTGeometry3D.X = Geometry3D.X(1:(Geometry3D.Ny-1),1:(Geometry3D.Nx-1),1:(Geometry3D.Nz-1)); % Restriction on the spatial grid on the x direction
FFTGeometry3D.Y = Geometry3D.Y(1:(Geometry3D.Ny-1),1:(Geometry3D.Nx-1),1:(Geometry3D.Nz-1)); % Restriction on the spatial grid on the y direction
FFTGeometry3D.Z = Geometry3D.Z(1:(Geometry3D.Ny-1),1:(Geometry3D.Nx-1),1:(Geometry3D.Nz-1)); % Restriction on the spatial grid on the z direction
FFTGeometry3D.Nx = Geometry3D.Nx - 1; % Number of points in the x direction
FFTGeometry3D.Ny = Geometry3D.Ny - 1; % Number of points in the y direction
FFTGeometry3D.Nz = Geometry3D.Nz - 1; % Number of points in the z direction
FFTGeometry3D.N3 = FFTGeometry3D.Nx*FFTGeometry3D.Ny*FFTGeometry3D.Nz; % Total number of interior points