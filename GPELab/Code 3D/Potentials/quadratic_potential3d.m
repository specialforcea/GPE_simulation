%% Computation of a quadratic potential
%% INPUTS:
%%          gamma_x,gamma_y,gamma_z: Potential parameters (double)
%%          X,Y,Z: Computational grids (matrix)
%% OUTPUT:
%%          Potential : values of the potential over the computational grid (matrix)

function [Potential] = quadratic_potential3d(gamma_x, gamma_y,gamma_z,X,Y,Z)
Potential = (gamma_x^2*X.^2 + gamma_y^2*Y.^2+ gamma_z^2*Z.^2)/2; % Computing the quadratic potential