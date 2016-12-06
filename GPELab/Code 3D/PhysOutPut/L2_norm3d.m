%% Computation of the L2_norm
%% INPUTS:
%%          phi: Wave function whose L2 norm is computed (matrix)
%%          Geometry3D: Structure containing variables concerning the geometry of the problem in 3D (structure) (see Geometry3D_Var3d.m)
%% OUTPUT:
%%          L2norm: L2 norm of the function phi (double)

function [L2norm]= L2_norm3d(phi,Geometry3D)
L2norm = sqrt(Geometry3D.dx*Geometry3D.dy*Geometry3D.dz)*sqrt(sum(sum(sum(abs(phi).^2)))); % Computation of the L2 norm