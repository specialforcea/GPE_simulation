%% Computation of the gradient in the x direction with a four points finite difference scheme and Dirichlet conditions
%% INPUT:
%%          FDGeometry2D: Structure containing variables concerning the geometry of the problem in 2D in the FD context (structure) (see FDGeometry2D_Var2d.m)
%% OUTPUT:
%%          Gx: Discretized gradient in the x direction with the FD scheme (matrix)

function [Gx] = Four_Points_Gradientx2d(FDGeometry2D)
% Diagonals for the derivative in the x direction
h_x = 1/(12*FDGeometry2D.dx);
Diag_x1 = -8*h_x*spdiags(ones(FDGeometry2D.N2,1),FDGeometry2D.Ny,FDGeometry2D.N2,FDGeometry2D.N2); % Computing the matrix corresponding to the first half of the gradient in the x direction
Diag_x2 = h_x*spdiags(ones(FDGeometry2D.N2,1),2*FDGeometry2D.Ny,FDGeometry2D.N2,FDGeometry2D.N2); % Computing the matrix corresponding to the second half of the gradient in the x direction
Gx = Diag_x1 + Diag_x2 - Diag_x1' - Diag_x2'; % Assembling the gradient operator