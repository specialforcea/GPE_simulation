%% Computation of the root mean square of a two functions
%% INPUTS:
%%          phi,alpha: functions whose mean square is computed (matrix)
%%          Geometry3D: Structure containing variables concerning the geometry of the problem in 3D (structure) (see Geometry3D_Var3d.m)
%% OUTPUT:
%%          alpha_rms: Alpha root mean square of phi (double)
%% FUNCTIONS USED:
%%          L2_norm3d: To compute the root mean square (line 11)

function [alpha_rms]= alpha_rms3d(phi,alpha,Geometry3D)
alpha_rms = L2_norm3d(alpha.*phi,Geometry3D); % Computation of the root mean square of alpha and phi