%% Computation of a quadruple well trapping potential
%% INPUT:
%%          gamma_x,gamma_y,gamma_z,alpha_x,alpha_y,aplha_z: potential parameters (double)
%%          X,Y,Z: Computational grids (matrix)
%% OUTPUT:
%%          Potential : values of the potential over the computational grid (matrix)

function [Potential] = quadruple_well_trapping_potential3d(gamma_x,gamma_y,gamma_z,alpha_x,alpha_y,alpha_z,X,Y,Z)
Potential= alpha_x*X.^4 + alpha_y*Y.^4 + alpha_z*Z.^4 + gamma_x*X.^2 + gamma_y*Y.^2 + gamma_z*Z.^2;  % Computing the quadruple well trapping potential