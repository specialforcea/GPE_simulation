%% Computation of the Laplacian with a five points finite difference scheme with Dirichlet conditions
%% INPUT:
%%          FDGeometry2D: Structure containing variables concerning the geometry of the problem in 2D in the FD context (structure) (see FDGeometry2D_Var2d.m)
%% OUTPUT:
%%          A: Discretized Laplacian with the FD scheme (matrix)

function [A] = Five_Points_Laplacian2d(FDGeometry2D)
%% Computing the diagonals
h_central = -2*(FDGeometry2D.dx^2+FDGeometry2D.dy^2)/(FDGeometry2D.dx*FDGeometry2D.dy)^2;
Diag_central = spdiags([ones(FDGeometry2D.N2,1)],0,FDGeometry2D.N2,FDGeometry2D.N2)*h_central; %Computing the central diagonal matrice of the discretized Laplacian
h_y = 1/FDGeometry2D.dy^2;
Diag_y = repmat([0;ones(FDGeometry2D.Ny-1,1)],FDGeometry2D.Nx,1)*h_y; %Computing the vector corresponding to half of the gradient in the y direction
h_x = 1/FDGeometry2D.dx^2;
Diag_x = ones(FDGeometry2D.N2,1)*h_x;%Computing the vector corresponding to half of the gradient in the x direction
Diag_xy = spdiags(Diag_y,1,FDGeometry2D.N2,FDGeometry2D.N2)+spdiags(Diag_x,FDGeometry2D.Ny,FDGeometry2D.N2,FDGeometry2D.N2); %Computing the matrice with the two half gradients
%% Computing the discretized Laplacian
A = Diag_central+Diag_xy+Diag_xy'; %Adding the two gradients