%% Creation of the 2D geometry structure for the FFT
%% INPUT:
%%          Geometry2D: Structure containing variables concerning the 2D geometry (structure) (see Geometry2D_Var2d.m)
%% OUTPUT:
%%          FFTGeometry2D: Structure containing variables concerning the 2D geometry for the FFT (structure)

function [FFTGeometry2D] = FFTGeometry2D_Var2d(Geometry2D)
%% Initialization
FFTGeometry2D = Geometry2D; % Copying the 2D geometry

%% Preparation for the FFT of the 2D geometry
FFTGeometry2D.X = Geometry2D.X(1:(Geometry2D.Ny-1),1:(Geometry2D.Nx-1)); % Restriction on the spatial grid on the x direction
FFTGeometry2D.Y = Geometry2D.Y(1:(Geometry2D.Ny-1),1:(Geometry2D.Nx-1)); % Restriction on the spatial grid on the y direction
FFTGeometry2D.Nx = Geometry2D.Nx - 1; % Number of points in the x direction
FFTGeometry2D.Ny = Geometry2D.Ny - 1; % Number of points in the y direction
FFTGeometry2D.N2 = FFTGeometry2D.Nx*FFTGeometry2D.Ny; % Total number of interior points