%% Computation of the Thomas-Fermi approximation
%% INPUTS:
%%          gamma_x: Potential parameters (double)
%%          Beta: Coefficient in front of the nonlinearity (double)
%%          Potential: Potential of the problem (matrix)
%% OUTPUT:
%%          phi_0: Thomas-Fermi approximation (matrix)

function [phi_0] = Thomas_Fermi1d(gamma_x,Beta,Potential)
phi_0 = real(sqrt((1/2*sqrt((3*Beta*gamma_x/2)) - Potential)/Beta)); % Computing the Thomas-Fermi approximation