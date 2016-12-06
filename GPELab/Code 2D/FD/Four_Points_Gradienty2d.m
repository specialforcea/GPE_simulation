%% Computation of the gradient in the y direction with a four points finite difference scheme and Dirichlet conditions
%% INPUT:
%%          FDGeometry2D: Structure containing variables concerning the geometry of the problem in 2D in the FD context (structure) (see FDGeometry2D_Var2d.m)
%% OUTPUT:
%%          Gy: Discretized gradient in the y direction with the FD scheme (matrix)

function [Gy] = Four_Points_Gradienty2d(FDGeometry2D);
% Diagonals for the derivative in the y direction
h_y = 1/(12*FDGeometry2D.dy);
Diag_y1 = -8*h_y*spdiags(repmat([0;ones(FDGeometry2D.Ny-1,1)],FDGeometry2D.Nx,1), 1 ,FDGeometry2D.N2,FDGeometry2D.N2); %Computing the vector corresponding to the first half of the gradient in the y direction
Diag_y2 = h_y*spdiags(repmat([0;0;ones(FDGeometry2D.Ny-2,1)],FDGeometry2D.Nx,1), 2 ,FDGeometry2D.N2,FDGeometry2D.N2); %Computing the vector corresponding to the second half of the gradient in the y direction
Gy = Diag_y1 + Diag_y2 - Diag_y1' - Diag_y2'; %Assembling the gradient operator