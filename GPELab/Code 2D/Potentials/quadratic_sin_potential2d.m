%% Computation of a harmonic plus sinusoidal terms potential
%% INPUTS:
%%          gamma_x,gamma_y: Potential parameters (double)
%%          a: Intensities of the sinusoidal terms (vector)
%%          d: Frequency parameter (vector)
%%          X,Y: Computational grids (matrix)
%% OUTPUT:
%%          Potential : values of the potential over the computational grid (matrix)

function [Potential] = quadratic_sin_potential2d(gamma_x, gamma_y,a,d,X,Y)
Potential = (gamma_x*X.^2 + gamma_y*Y.^2 + a(1).*sin(pi.*X/d(1)).^2 + a(2)*sin(pi.*Y/d(2)).^2)/2; % Computing the harmonic plus sinusoidal terms potential