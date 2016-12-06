%% Computation of the angular momentum using the Fourier transform
%% INPUTS:
%%          phi: Wave function (matrix)
%%          FFTGeometry2D: Structure containing variables concerning the geometry of the problem in 2D in the FFT context (structure) (see FFTGeometry2D_Var2d.m)
%%          FFTOperators2D: Structure containing the derivative FFT operators (structure) (see FFTOperators2D_Var2d.m)
%% OUTPUT:
%%          Angular_momentum: Angular momentum of the wave function (double)

function [Angular_momentum] = Angular_momentum_Fourier2d(phi, FFTGeometry2D, FFTOperators2D)

%% Derivation using the Fourier transform
Gradx_phi = ifft(FFTOperators2D.Gx.*fft(phi,[],2),[],2); % First order derivating the component's wave function on the x direction via FFT and iFFT
Grady_phi = ifft(FFTOperators2D.Gy.*fft(phi,[],1),[],1); % First order derivating the component's wave function on the y direction via FFT and iFFT

%% Computation of the local momentum and integration over space
Local_momentum = 1i*conj(phi).*(FFTGeometry2D.X.*Grady_phi - FFTGeometry2D.Y.*Gradx_phi); %Computation of the local momentum 
Angular_momentum = FFTGeometry2D.dx*FFTGeometry2D.dy*sum(sum(sum(Local_momentum))); %Integration over the space