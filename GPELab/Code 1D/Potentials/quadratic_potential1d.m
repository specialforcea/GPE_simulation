%% Computation of a quadratic potential
%% INPUTS:
%%          gamma_x: Potential parameters (double)
%%          X: Computational grids (vector)
%% OUTPUT:
%%          Potential : values of the potential over the computational grid (matrix)

function [Potential] = quadratic_potential1d(gamma_x,X)
Potential = gamma_x^2*X.^2/2; % Computing the quadratic potential