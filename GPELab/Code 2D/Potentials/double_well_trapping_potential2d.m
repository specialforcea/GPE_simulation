%% Computation of a quadruple well trapping potential
%% INPUTS:
%%          gamma_x,gamma_y,V_0,d: potential parameters (double)
%%          X,Y: Computational grids (matrix)
%% OUTPUT:
%%          Potential : values of the potential over the computational grid (matrix)

function [Potential] = double_well_trapping_potential2d(gamma_x, gamma_y,V_0,d,X,Y)
Potential= (1/2)*(gamma_x*X.^2 + gamma_y*Y.^2) + V_0*exp(-X.^2/(2*d));  % Computing the quadruple well trapping potential