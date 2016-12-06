%% Computation of a quadruple well trapping potential
%% INPUTS:
%%          gamma_x,alpha_x: potential parameters (double)
%%          X: Computational grid (vecotr)
%% OUTPUT:
%%          Potential : values of the potential over the computational grid (vector)

function [Potential] = quadruple_well_trapping_potential1d(gamma_x,alpha_x,X)
Potential= alpha_x*X.^4+ gamma_x*X.^2;  % Computing the quadruple well trapping potential