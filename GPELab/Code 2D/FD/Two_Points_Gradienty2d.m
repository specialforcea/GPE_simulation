%% Computation of the gradient in the y direction with a two points finite difference scheme and Dirichlet conditions
%% INPUT:
%%          FDGeometry2D: Structure containing variables concerning the geometry of the problem in 2D in the FD context (structure) (see FDGeometry2D_Var2d.m)
%% OUTPUT:
%%          Gy: Discretized gradient in the y direction with the FD scheme (matrix)

function [Gy] = Two_Points_Gradienty2d(FDGeometry2D);
N2=FDGeometry2D.Nx*FDGeometry2D.Ny; %Number of internal points of the 2D spatial grid
Diagy=repmat([0;ones(FDGeometry2D.Ny-1,1)],FDGeometry2D.Nx,1); %Computing the vector corresponding to the diagonal for the gradient operator
Diagy=spdiags(Diagy,1,N2,N2); %Computing the diagonal for the gradient operator
Gy=(1/(2*FDGeometry2D.dy))*(Diagy-Diagy'); %Assembling the gradient operator