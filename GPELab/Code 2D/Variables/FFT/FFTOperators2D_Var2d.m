%% Creation of the 2D operators structure for the FFT
%% INPUT:
%%          FFTGeometry2D: Structure containing variables concerning the 2D geometry for the FFT (structure) (see FFTGeometry2D_Var2d.m)
%% OUTPUT:
%%          FFTOperators2D: Structure containing the derivative FFT operators (structure)
%% FUNCTIONS USED:
%%          Delta_Fourier2d: To compute the second derivative FFT operators (line 13)
%%          Grad_Fourier2d: To compute the first derivative FFT operators (line 14)

function [FFTOperators2D] = FFTOperators2D_Var2d(FFTGeometry2D)

%% Computation of the FFT derivative operators
[FFTOperators2D.Dx,FFTOperators2D.Dy] = Delta_Fourier2d(FFTGeometry2D); % Computing the second derivative FFT operators
[FFTOperators2D.Gx,FFTOperators2D.Gy] = Grad_Fourier2d(FFTGeometry2D); % Computing the first derivative FFT operators