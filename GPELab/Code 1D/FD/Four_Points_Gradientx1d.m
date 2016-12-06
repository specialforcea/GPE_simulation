%% Computation of the gradient in the x direction with a four points finite difference scheme and Dirichlet conditions
%% INPUT:
%%          FDGeometry1D: Structure containing variables concerning the geometry of the problem in 1D in the FD context (structure) (see FDGeometry1D_Var1d.m)
%% OUTPUT:
%%          Gx: Discretized gradient in the x direction with the FD scheme (matrix)

function [Gx] = Four_Points_Gradientx1d(FDGeometry1D)
%% Computing the distretized gradient
h = 1/(12*FDGeometry1D.dx);
Diag_x1 = -8*h*spdiags(ones(FDGeometry1D.Nx,1),1,FDGeometry1D.Nx,FDGeometry1D.Nx); % Computing the first diagonal of the gradient operator
Diag_x2 = h*spdiags(ones(FDGeometry1D.Nx,1),2,FDGeometry1D.Nx,FDGeometry1D.Nx); % Computing the second diagonal of the gradient operator
Gx = Diag_x1 - Diag_x1' + Diag_x2 - Diag_x2'; % Assembling the gradient operator