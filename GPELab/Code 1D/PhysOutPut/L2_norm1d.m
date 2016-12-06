%% Computation of the L2_norm
%% INPUTS:
%%          phi: Wave function whose L2 norm is computed (matrix)
%%          Geometry1D: Structure containing variables concerning the geometry of the problem in 1D (structure) (see Geometry1D_Var1d.m)
%% OUTPUT:
%%          L2norm: L2 norm of the function phi (double)

function [L2norm]= L2_norm1d(phi,Geometry1D)
L2norm = sqrt(Geometry1D.dx)*sqrt(sum(abs(phi).^2)); % Computation of the L2 norm