%% Computation of a quadratic plus quartic potential 
%% INPUTS:
%%          gamma_x,alpha,kappa: Potential parameters (double)
%%          X: Computational grids (vector)
%% OUTPUT:
%%          Potential : values of the potential over the computational grid (matrix)

function [Potential] = quadratic_plus_quartic_potential1d(gamma_x,alpha,kappa,X)
Potential = (1-alpha)*gamma_x*X.^2/2 + kappa*(gamma_x*X.^2 /2).^2; % Computing the quadratic plus quartic potential