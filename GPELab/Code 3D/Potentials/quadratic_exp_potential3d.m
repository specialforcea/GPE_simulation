%% Computation of a harmonic plus exponential term potential (far-blue detuned Gaussian laser beam)
%% INPUTS:
%%          gamma_x,gamma_y,gamma_z: Quadratic potential parameters (double)
%%          X0,Y0,Z0: Coordinates of the center of the exponential potential (double)
%%          d,w_0: Exponential potential parameters (double)
%%          X,Y,Z: Computational grids (matrix)
%% OUTPUT:
%%          Potential : values of the potential over the computational grid (matrix)

function [Potential] = quadratic_exp_potential3d(gamma_x, gamma_y, gamma_z, X0, Y0, Z0, d, w_0,X,Y,Z)
Potential = (gamma_x^2*X.^2 + gamma_y^2*Y.^2 + gamma_z^2*Z.^2)/2 + w_0.*exp(-((X-X0).^2+(Y-Y0).^2+(Z-Z0).^2)/d^2); % Computing the harmonic plus exponential term potential