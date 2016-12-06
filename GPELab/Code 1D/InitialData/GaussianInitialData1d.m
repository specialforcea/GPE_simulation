%% Computation of the centered gaussian
%% INPUTS:
%%          Geometry1D: Structure containing variables concerning the geometry of the problem in 1D (structure) (see Geometry1D_Var1d.m)
%%          Physics1D: Structure containing variables concerning the physics of the problem in 1D (structure) (see Physics1D_Var1d.m)
%%          gamma_x: Potential parameters (double)
%%          X0: Coordinates of the center of the gaussian (double)
%% OUTPUT:
%%          phi_0: Centered gaussian (matrix)

function [phi_0] = GaussianInitialData1d(Geometry1D, gamma_x,  X0)
phi_0 = (gamma_x)^(1/2).*exp(-gamma_x*(Geometry1D.X-X0).^2/2)/sqrt(pi); % Computing the centered gaussian
phi_0 = phi_0/L2_norm1d(phi_0,Geometry1D); % Normalizing the centered gaussian