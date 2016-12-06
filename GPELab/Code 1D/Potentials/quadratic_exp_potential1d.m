%% Computation of a harmonic plus exponential term potential (far-blue detuned Gaussian laser beam)
%% INPUTS:
%%          gamma_x: Quadratic potential parameters (double)
%%          X0: Coordinates of the center of the exponential potential (double)
%%          d,w_0: Exponential potential parameters (double)
%%          X: Computational grids (vector)
%% OUTPUT:
%%          Potential : values of the potential over the computational grid (matrix)

function [Potential] = quadratic_exp_potential1d(gamma_x,X0, d, w_0,X)
Potential = gamma_x^2*X.^2 /2 + w_0.*exp(-(X-X0).^2/d^2); % Computing the harmonic plus exponential term potential