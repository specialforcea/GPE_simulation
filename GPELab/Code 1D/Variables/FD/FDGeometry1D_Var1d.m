%% Creation of the 1D geometry structure for the FD
%% INPUT:
%%          Geometry1D: Structure containing variables concerning the 1D geometry (structure) (see Geometry1D_Var1d.m)
%% OUTPUT:
%%          FDGeometry1D: Structure containing variables concerning the 1D geometry for the FD (structure)

function [FDGeometry1D] = FDGeometry1D_Var1d(Geometry1D)
%% Initialization
FDGeometry1D = Geometry1D; % Copying the 1D geometry

%% Preparation for the FD of the 1D geometry
FDGeometry1D.X = Geometry1D.X(2:(Geometry1D.Nx-1)); %Restriction on the spatial grid on the x direction
FDGeometry1D.Nx = Geometry1D.Nx - 2; %Number of interior points in the x direction