%% Computation of a harmonic plus exponential term potential (far-blue detuned Gaussian laser beam)
%% INPUTS:
%%          gamma_x,gamma_y: Quadratic potential parameters (double)
%%          X0,Y0: Coordinates of the center of the exponential potential (double)
%%          d,w_0: Exponential potential parameters (double)
%%          X,Y: Computational grids (matrix)
%% OUTPUT:
%%          Potential : values of the potential over the computational grid (matrix)

function [Potential] = quadratic_exp_potential2d(gamma_x, gamma_y,X0, Y0, d, w_0,X,Y)
Potential = (gamma_x^2*X.^2 + gamma_y^2*Y.^2)/2 + w_0.*exp(-((X-X0).^2+(Y-Y0).^2)/d^2); % Computing the harmonic plus exponential term potential