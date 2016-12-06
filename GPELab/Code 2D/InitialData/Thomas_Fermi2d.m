%% Computation of the Thomas-Fermi approximation
%% INPUTS:
%%          gamma_x,gamma_y: Potential parameters (double)
%%          Beta: Coefficient in front of the nonlinearity (double)
%%          Potential: Potential of the problem (matrix)
%% OUTPUT:
%%          phi_0: Thomas-Fermi approximation (matrix)

function [phi_0] = Thomas_Fermi2d(gamma_x,gamma_y,Beta,Potential)
phi_0 = real(sqrt((1/2*sqrt((4*Beta*gamma_x*gamma_y/pi)) - Potential)/Beta)); % Computing the Thomas-Fermi approximation