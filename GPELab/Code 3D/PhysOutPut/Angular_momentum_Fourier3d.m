%% Computation of the angular momentum using the Fourier transform
%% INPUTS:
%%          phi: Wave function (matrix)
%%          FFTGeometry3D: Structure containing variables concerning the geometry of the problem in 3D in the FFT context (structure) (see FFTGeometry3D_Var3d.m)
%%          FFTOperators3D: Structure containing the derivative FFT operators (structure) (see FFTOperators3D_Var3d.m)
%% OUTPUT:
%%          Angular_momentum: Angular momentum of the wave function (double)

function [Angular_momentum] = Angular_momentum_Fourier3d(phi, FFTGeometry3D, FFTOperators3D)

%% Derivation using the Fourier transform
Gradx_phi = ifft(FFTOperators3D.Gx.*fft(phi,[],2),[],2); % First order derivating the component's wave function on the x direction via FFT and iFFT
Grady_phi = ifft(FFTOperators3D.Gy.*fft(phi,[],1),[],1); % First order derivating the component's wave function on the y direction via FFT and iFFT
Gradz_phi = ifft(FFTOperators3D.Gz.*fft(phi,[],3),[],3); % First order derivating the component's wave function on the z direction via FFT and iFFT

%% Computation of the local momentum and integration over space
Local_momentum = 1i*conj(phi).*((FFTGeometry3D.X-FFTGeometry3D.Z).*Grady_phi + (FFTGeometry3D.Z - FFTGeometry3D.Y).*Gradx_phi + (FFTGeometry3D.Z + FFTGeometry3D.X).*Gradz_phi); %Computation of the local momentum 
Angular_momentum = FFTGeometry3D.dx*FFTGeometry3D.dy*FFTGeometry3D.dz*sum(sum(sum(Local_momentum))); %Integration over the space