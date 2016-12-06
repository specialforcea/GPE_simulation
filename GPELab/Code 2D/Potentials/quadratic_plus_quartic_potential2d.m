%% Computation of a quadratic plus quartic potential 
%% INPUTS:
%%          gamma_x,gamma_y,alpha,kappa: Potential parameters (double)
%%          X,Y: Computational grids (matrix)
%% OUTPUT:
%%          Potential : values of the potential over the computational grid (matrix)

function [Potential] = quadratic_plus_quartic_potential2d(gamma_x, gamma_y,alpha,kappa,X,Y)
Potential = (1-alpha)*(gamma_x*X.^2 + gamma_y*Y.^2)/2 + kappa*((gamma_x*X.^2 + gamma_y*Y.^2)/2).^2; % Computing the quadratic plus quartic potential