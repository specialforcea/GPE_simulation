%% Computation of the Thomas-Fermi approximation
%% INPUTS:
%%          gamma_x,gamma_y,gamma_z: Potential parameters (double)
%%          Beta: Coefficient in front of the nonlinearity (double)
%%          Potential: Potential of the problem (matrix)
%% OUTPUT:
%%          phi_0: Thomas-Fermi approximation (matrix)

function [phi_0] = Thomas_Fermi3d(gamma_x,gamma_y,gamma_z,Beta,Potential)
phi_0 = real(sqrt((1/2*(15*Beta*gamma_x*gamma_y*gamma_z/4/pi)^(2/5) - Potential)/Beta)); % Computing the Thomas-Fermi approximation