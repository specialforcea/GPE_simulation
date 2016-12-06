%% Computation of the gradient in the x direction with a two points finite difference scheme and Dirichlet conditions
%% INPUT:
%%          FDGeometry1D: Structure containing variables concerning the geometry of the problem in 1D in the FD context (structure) (see FDGeometry1D_Var1d.m)
%% OUTPUT:
%%          Gx: Discretized gradient in the x direction with the FD scheme (matrix)

function [Gx] = Two_Points_Gradientx1d(FDGeometry1D)
%% Computing the distretized gradient
Diagx=spdiags([ones(FDGeometry1D.Nx,1)],1,FDGeometry1D.Nx,FDGeometry1D.Nx); %Computing the diagonal of the gradient operator
Gx=(1/(2*FDGeometry1D.dx))*(Diagx-Diagx'); %Assembling the gradient operator