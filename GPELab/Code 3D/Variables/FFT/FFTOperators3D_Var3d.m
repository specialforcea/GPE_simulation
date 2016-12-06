%% Creation of the 3D operators structure for the FFT
%% INPUT:
%%          FFTGeometry3D: Structure containing variables concerning the 3D geometry for the FFT (structure) (see FFTGeometry3D_Var3d.m)
%% OUTPUT:
%%          FFTOperators3D: Structure containing the derivative FFT operators (structure)
%% FUNCTIONS USED:
%%          Delta_Fourier3d: To compute the second derivative FFT operators (line 13)
%%          Grad_Fourier3d: To compute the first derivative FFT operators (line 14)

function [FFTOperators3D] = FFTOperators3D_Var3d(FFTGeometry3D)

%% Computation of the FFT derivative operators
[FFTOperators3D.Dx,FFTOperators3D.Dy,FFTOperators3D.Dz]=Delta_Fourier3d(FFTGeometry3D); % Computing the second derivative FFT operators
[FFTOperators3D.Gx,FFTOperators3D.Gy,FFTOperators3D.Gz]=Grad_Fourier3d(FFTGeometry3D); % Computing the first derivative FFT operators