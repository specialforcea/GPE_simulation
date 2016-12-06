%% Computation of a quadruple well trapping potential
%% INPUTS:
%%          gamma_x,gamma_y,alpha_x,alpha_y: potential parameters (double)
%%          X,Y: Computational grids (matrix)
%% OUTPUT:
%%          Potential : values of the potential over the computational grid (matrix)

function [Potential] = quadruple_well_trapping_potential2d(gamma_x, gamma_y,alpha_x,alpha_y,X,Y)
Potential= alpha_x*X.^4 + alpha_y*Y.^4 + gamma_x*X.^2 + gamma_y*Y.^2;  % Computing the quadruple well trapping potential