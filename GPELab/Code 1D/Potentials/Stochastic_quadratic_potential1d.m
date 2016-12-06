%% Computation of a stochastic quadratic potential
%% INPUTS:
%%          gamma_x: Potential parameters (double)
%%          X: Computational grids (vector)
%% OUTPUT:
%%          Potential : values of the potential over the computational grid (matrix)

function [Potential] = Stochastic_quadratic_potential1d(Method,gamma_x,X)
Potential = sqrt(Method.Deltat)^(-1)*randn*gamma_x^2*X.^2/2; % Computing the quadratic potential