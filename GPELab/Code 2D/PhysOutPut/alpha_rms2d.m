%% Computation of the root mean square of a two functions
%% INPUTS:
%%          phi,alpha: functions whose mean square is computed (matrix)
%%          Geometry2D: Structure containing variables concerning the geometry of the problem in 2D (structure) (see Geometry2D_Var2d.m)
%% OUTPUT:
%%          alpha_rms: Alpha root mean square of phi (double)
%% FUNCTIONS USED:
%%          L2_norm2d: To compute the root mean square (line 11)

function [alpha_rms]= alpha_rms2d(phi,alpha,Geometry2D)
alpha_rms = L2_norm2d(alpha.*phi,Geometry2D); % Computation of the root mean square of alpha and phi