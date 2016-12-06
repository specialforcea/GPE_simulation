%% Creation of the 1D geometry structure for the FFT
%% INPUT:
%%          Geometry1D: Structure containing variables concerning the 1D geometry (structure) (see Geometry1D_Var1d.m)
%% OUTPUT:
%%          FFTGeometry1D: Structure containing variables concerning the 1D geometry for the FFT (structure)

function [FFTGeometry1D] = FFTGeometry1D_Var1d(Geometry1D)
%% Initialization
FFTGeometry1D = Geometry1D; % Copying the 1D geometry

%% Preparation for the FFT of the 1D geometry
FFTGeometry1D.X = Geometry1D.X(1:(Geometry1D.Nx-1)); % Restriction on the spatial grid on the x direction
FFTGeometry1D.Nx = Geometry1D.Nx - 1; % Number of points in the x direction