%% Computation of the centered gaussian
%% INPUTS:
%%          Geometry3D: Structure containing variables concerning the geometry of the problem in 3D (structure) (see Geometry3D_Var3d.m)
%%          Physics3D: Structure containing variables concerning the physics of the problem in 3D (structure) (see Physics3D_Var3d.m)
%%          gamma_x,gamma_y,gamma_z: Potential parameters (double)
%%          X0,Y0,Z0: Coordinates of the center of the gaussian (double)
%% OUTPUT:
%%          phi_0: Centered gaussian (matrix)

function [phi_0] = GaussianInitialData3d(Geometry3D, Physics3D, gamma_x, gamma_y, gamma_z,  X0, Y0, Z0)
phi_0 = (gamma_x*gamma_y*gamma_z)^(1/4)*((1-sum(Physics3D.Omega)).*exp(-(gamma_x*(Geometry3D.X-X0).^2+gamma_y*(Geometry3D.Y-Y0).^2+gamma_z*(Geometry3D.Z-Z0).^2)/2)/sqrt(pi)) ...
        + Physics3D.Omega(3).*((gamma_x*(Geometry3D.X-X0)+1i*gamma_y*(Geometry3D.Y-Y0)).*exp(-(gamma_x*(Geometry3D.X-X0).^2+gamma_y*(Geometry3D.Y-Y0).^2+gamma_z*(Geometry3D.Z-Z0).^2)/2)/sqrt(pi))...
        + Physics3D.Omega(1).*((gamma_y*(Geometry3D.Y-Y0)+1i*gamma_z*(Geometry3D.Z-Z0)).*exp(-(gamma_x*(Geometry3D.X-X0).^2+gamma_y*(Geometry3D.Y-Y0).^2+gamma_z*(Geometry3D.Z-Z0).^2)/2)/sqrt(pi))...
        + Physics3D.Omega(2).*((gamma_z*(Geometry3D.Z-Z0)+1i*gamma_x*(Geometry3D.X-X0)).*exp(-(gamma_x*(Geometry3D.X-X0).^2+gamma_y*(Geometry3D.Y-Y0).^2+gamma_z*(Geometry3D.Z-Z0).^2)/2)/sqrt(pi)); % Computing the centered gaussian
phi_0 = phi_0./L2_norm3d(phi_0,Geometry3D); % Normalizing the centered gaussian