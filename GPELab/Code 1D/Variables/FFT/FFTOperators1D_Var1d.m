%% Creation of the 1D operators structure for the FFT
%% INPUT:
%%          FFTGeometry1D: Structure containing variables concerning the 1D geometry for the FFT (structure) (see FFTGeometry1D_Var1d.m)
%% OUTPUT:
%%          FFTOperators1D: Structure containing the derivative FFT operators (structure)
%% FUNCTIONS USED:
%%          Delta_Fourier1d: To compute the second derivative FFT operators (line 13)
%%          Grad_Fourier1d: To compute the first derivative FFT operators (line 14)

function [FFTOperators1D] = FFTOperators1D_Var1d(FFTGeometry1D)

%% Computation of the FFT derivative operators
[FFTOperators1D.Dx]=Delta_Fourier1d(FFTGeometry1D); % Computing the second derivative FFT operators
[FFTOperators1D.Gx]=Grad_Fourier1d(FFTGeometry1D); % Computing the first derivative FFT operators