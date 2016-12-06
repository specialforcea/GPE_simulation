%% Computation of the Laplacian with a three points finite difference scheme with Dirichlet conditions
%% INPUT:
%%          FDGeometry1D: Structure containing variables concerning the geometry of the problem in 1D in the FD context (structure) (see FDGeometry1D_Var1d.m)
%% OUTPUT:
%%          A: Discretized Laplacian with the FD scheme (matrix)

function [A] = Three_Points_Laplacian1d(FDGeometry1D);
%% Computing the diagonals
h_central = -2/(FDGeometry1D.dx)^2;
Diag_central = spdiags([ones(FDGeometry1D.Nx,1)],0,FDGeometry1D.Nx,FDGeometry1D.Nx)*h_central; % Computing the central diagonal matrice of the discretized Laplacian
h_x = 1/(FDGeometry1D.dx)^2;
Diag_x = ones(FDGeometry1D.Nx,1)*h_x; % Computing the vector corresponding to half of the gradient in the x direction
Diag_x = spdiags(Diag_x,1,FDGeometry1D.Nx,FDGeometry1D.Nx); % Computing the matrice with the half gradient in the x direction
%% Computing the discretized Laplacian
A = Diag_central+Diag_x+Diag_x'; %Adding the two gradients