%% Computation of the Laplacian with a ten points finite difference scheme with Dirichlet conditions
%% INPUT:
%%          FDGeometry2D: Structure containing variables concerning the geometry of the problem in 2D in the FD context (structure) (see FDGeometry2D_Var2d.m)
%% OUTPUT:
%%          A: Discretized Laplacian with the FD scheme (matrix)

function [A] = Ten_Points_Laplacian2d(FDGeometry2D)
%% Computing the diagonals
h_central = -30/(12*FDGeometry2D.dx^2) - 30/(12*FDGeometry2D.dy^2);
Diag_central = spdiags([ones(FDGeometry2D.N2,1)],0,FDGeometry2D.N2,FDGeometry2D.N2)*h_central; %Computing the central diagonal matrice of the discretized Laplacian
% Diagonals for the derivative in the x direction
h_x = 1/(12*FDGeometry2D.dx^2);
Diag_x1 = 16*h_x*spdiags(ones(FDGeometry2D.N2,1),FDGeometry2D.Ny,FDGeometry2D.N2,FDGeometry2D.N2); % Computing the matrix corresponding to the first half of the gradient in the x direction
Diag_x2 = -h_x*spdiags(ones(FDGeometry2D.N2,1),2*FDGeometry2D.Ny,FDGeometry2D.N2,FDGeometry2D.N2); % Computing the matrix corresponding to the second half of the gradient in the x direction
% Diagonals for the derivative in the y direction
h_y = 1/(12*FDGeometry2D.dy^2);
Diag_y1 = 16*h_y*spdiags(repmat([0;ones(FDGeometry2D.Ny-1,1)],FDGeometry2D.Nx,1), 1 ,FDGeometry2D.N2,FDGeometry2D.N2); %Computing the vector corresponding to the first half of the gradient in the y direction
Diag_y2 = -h_y*spdiags(repmat([0;0;ones(FDGeometry2D.Ny-2,1)],FDGeometry2D.Nx,1), 2 ,FDGeometry2D.N2,FDGeometry2D.N2); %Computing the vector corresponding to the second half of the gradient in the y direction
%% Computing the discretized Laplacian
A = Diag_central + Diag_x1 + Diag_x1' + Diag_x2 + Diag_x2' + Diag_y1 + Diag_y1' + Diag_y2 + Diag_y2'; %Adding the two gradients