function [B] = Stationary_Gaussian_Field1d(Geometry1D,varargin)
%% Analysis of the inputs
Analyse_Var = inputParser; % Creating the parser
Analyse_Var.addOptional('Cov_fun', @(X) exp(-X.^2)); % Required input 'Cov_fun' corresponding to a covariance function

%% Computing the geometry for the FFT
FFTGeometry1D = FFTGeometry1D_Var1d(Geometry1D);

%% Parsing inputs and storing inputs
% Parsing inputs
Analyse_Var.parse(varargin{:}); % Analysing the inputs
Cov_fun = Analyse_Var.Results.Cov_fun; %Storing the 'Method' input

% Computing a random normal matrix
RandN = randn(FFTGeometry1D.Nx,1);

% Computing the gaussian field using the spectral representation
Cov = Cov_fun(FFTGeometry1D.X);
FFTRandN = fft(RandN);
FFTCov = fft(fftshift(Cov));
B = real(ifft(sqrt(FFTCov).*FFTRandN));