%% Computation of the L2_norm
%% INPUTS:
%%          phi: Wave function whose L2 norm is computed (matrix)
%%          Geometry2D: Structure containing variables concerning the geometry of the problem in 2D (structure) (see Geometry2D_Var2d.m)
%% OUTPUT:
%%          L2norm: L2 norm of the function phi (double)

function [L2norm]= L2_norm2d(phi,Geometry2D)
L2norm = sqrt(Geometry2D.dx*Geometry2D.dy)*sqrt(sum(sum(abs(phi).^2))); % Computation of the L2 norm