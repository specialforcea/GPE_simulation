%% Computation of the Laplacian with a five points finite difference scheme with Dirichlet conditions
%% INPUT:
%%          FDGeometry1D: Structure containing variables concerning the geometry of the problem in 1D in the FD context (structure) (see FDGeometry1D_Var1d.m)
%% OUTPUT:
%%          A: Discretized Laplacian with the FD scheme (matrix)

function [A] = Five_Points_Laplacian1d(FDGeometry1D)
%% Computing the diagonals
h = 1/(12*FDGeometry1D.dx^2);
Diag_central = -30*h*spdiags(ones(FDGeometry1D.Nx,1),0,FDGeometry1D.Nx,FDGeometry1D.Nx); %Computing the central diagonal matrice of the discretized Laplacian
Diag_x1 = 16*h*spdiags(ones(FDGeometry1D.Nx,1),1,FDGeometry1D.Nx,FDGeometry1D.Nx); %Computing the first diagonal matrice of the discretized Laplacian
Diag_x2 = -h*spdiags(ones(FDGeometry1D.Nx,1),2,FDGeometry1D.Nx,FDGeometry1D.Nx); %Computing the first diagonal matrice of the discretized Laplacian
%% Computing the discretized Laplacian
A = Diag_central + Diag_x1 + Diag_x1' + Diag_x2 + Diag_x2'; % Constructing the disctretized Laplacian