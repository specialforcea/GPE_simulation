%% Computation of the centered gaussian
%% INPUTS:
%%          Geometry2D: Structure containing variables concerning the geometry of the problem in 2D (structure) (see Geometry2D_Var2d.m)
%%          Physics2D: Structure containing variables concerning the physics of the problem in 2D (structure) (see Physics2D_Var2d.m)
%%          gamma_x,gamma_y: Potential parameters (double)
%%          X0,Y0: Coordinates of the center of the gaussian (double)
%% OUTPUT:
%%          phi_0: Centered gaussian (matrix)

function [phi_0] = GaussianInitialData2d(Geometry2D, Physics2D, gamma_x, gamma_y,  X0, Y0)
phi_0 = (gamma_x*gamma_y)^(1/4)*((1-Physics2D.Omega).*exp(-(gamma_x*(Geometry2D.X-X0).^2+gamma_y*(Geometry2D.Y-Y0).^2)/2)/sqrt(pi)) ...
        + Physics2D.Omega.*((gamma_x*(Geometry2D.X-X0)-1i*gamma_y*(Geometry2D.Y-Y0)).*exp(-(gamma_x*(Geometry2D.X-X0).^2+gamma_y*(Geometry2D.Y-Y0).^2)/2)/sqrt(pi)); % Computing the centered gaussian
phi_0 = phi_0./L2_norm2d(phi_0,Geometry2D); % Normalizing the centered gaussian