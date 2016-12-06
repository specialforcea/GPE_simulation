%% Computation of a harmonic plus sinusoidal terms potential
%% INPUTS:
%%          gamma_x Potential parameters (double)
%%          a: Intensities of the sinusoidal terms (double)
%%          d: Frequency parameter (double)
%%          X: Computational grids (vector)
%% OUTPUT:
%%          Potential : values of the potential over the computational grid (matrix)

function [Potential] = quadratic_sin_potential1d(gamma_x,a,d,X)
Potential = (gamma_x*X.^2 + a.*sin(pi.*X/d).^2)/2; % Computing the harmonic plus sinusoidal terms potential