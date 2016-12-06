%% Computation of a quadratic potential
%% INPUTS:
%%          gamma_x,gamma_y: Potential parameters (double)
%%          X,Y: Computational grids (matrix)
%% OUTPUT:
%%          Potential : values of the potential over the computational grid (matrix)

function [Potential] = quadratic_potential2d(gamma_x, gamma_y,X,Y)
Potential = (gamma_x^2*X.^2 + gamma_y^2*Y.^2)/2; % Computing the quadratic potential