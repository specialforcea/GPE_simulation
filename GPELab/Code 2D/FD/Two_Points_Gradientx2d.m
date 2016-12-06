%% Computation of the gradient in the x direction with a two points finite difference scheme and Dirichlet conditions
%% INPUT:
%%          FDGeometry2D: Structure containing variables concerning the geometry of the problem in 2D in the FD context (structure) (see FDGeometry2D_Var2d.m)
%% OUTPUT:
%%          Gx: Discretized gradient in the x direction with the FD scheme (matrix)

function [Gx] = Two_Points_Gradientx2d(FDGeometry2D)
N2=FDGeometry2D.Ny*FDGeometry2D.Nx; %Total number of internal grid points
Diagx=spdiags([ones(N2,1)],FDGeometry2D.Ny,N2,N2); %Computing the diagonal of the gradient operator
Gx=(1/(2*FDGeometry2D.dx))*(Diagx-Diagx'); %Assembling the gradient operator