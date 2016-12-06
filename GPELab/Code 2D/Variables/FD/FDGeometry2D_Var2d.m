%% Creation of the 2D geometry structure for the FD
%% INPUT:
%%          Geometry2D: Structure containing variables concerning the 2D geometry (structure) (see Geometry2D_Var2d.m)
%% OUTPUT:
%%          FDGeometry2D: Structure containing variables concerning the 2D geometry for the FD (structure)

function [FDGeometry2D] = FDGeometry2D_Var2d(Geometry2D)
%% Initialization
FDGeometry2D = Geometry2D; % Copying the 2D geometry

%% Preparation for the FD of the 2D geometry
FDGeometry2D.X = Geometry2D.X(2:(Geometry2D.Ny-1),2:(Geometry2D.Nx-1)); %Restriction on the spatial grid on the x direction
FDGeometry2D.Y = Geometry2D.Y(2:(Geometry2D.Ny-1),2:(Geometry2D.Nx-1)); %Restriction on the spatial grid on the y direction
FDGeometry2D.Nx = Geometry2D.Nx - 2; %Number of interior points in the x direction
FDGeometry2D.Ny = Geometry2D.Ny - 2; %Number of interior points in the y direction
FDGeometry2D.N2 = FDGeometry2D.Ny*FDGeometry2D.Nx; %Total number of interior points