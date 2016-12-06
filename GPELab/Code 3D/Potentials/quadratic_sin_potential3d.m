%% Computation of a harmonic plus sinusoidal terms potential
%% INPUTS:
%%          gamma_x,gamma_y,gamma_z: Potential parameters (double)
%%          a: Intensities of the sinusoidal terms (vector)
%%          d: Frequency parameter (vector)
%%          X,Y,Z: Computational grids (matrix)
%% OUTPUT:
%%          Potential : values of the potential over the computational grid (matrix)

function [Potential] = quadratic_sin_potential3d(gamma_x, gamma_y, gamma_z,a,d,X,Y,Z)
Potential = (gamma_x*X.^2 + gamma_y*Y.^2 + gamma_z*Z.^2 + a(1).*sin(pi.*X/d(1)).^2 + a(2)*sin(pi.*Y/d(2)).^2+ a(3).*sin(pi.*Z/d(3)).^2)/2; % Computing the harmonic plus sinusoidal terms potential