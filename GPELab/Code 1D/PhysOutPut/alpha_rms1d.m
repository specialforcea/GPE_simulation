%% Computation of the root mean square of a two functions
%% INPUTS:
%%          phi,alpha: functions whose mean square is computed (matrix)
%%          Geometry1D: Structure containing variables concerning the geometry of the problem in 1D (structure) (see Geometry1D_Var1d.m)
%% OUTPUT:
%%          alpha_rms: Alpha root mean square of phi (double)
%% FUNCTIONS USED:
%%          L2_norm1d: To compute the root mean square (line 11)

function [alpha_rms]= alpha_rms1d(phi,alpha,Geometry1D)
alpha_rms = L2_norm1d(alpha.*phi,Geometry1D); % Computation of the root mean square of alpha and phi