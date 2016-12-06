function [B] = Stationary_Gaussian_Field2d(Geometry2D,varargin)
%% Analysis of the inputs
Analyse_Var = inputParser; % Creating the parser
Analyse_Var.addOptional('Cov_fun', @(X,Y) exp(-X.^2-Y.^2)); % Required input 'Cov_fun' corresponding to a covariance function

%% Computing the geometry for the FFT
FFTGeometry2D = FFTGeometry2D_Var2d(Geometry2D);

%% Parsing inputs and storing inputs
% Parsing inputs
Analyse_Var.parse(varargin{:}); % Analysing the inputs
Cov_fun = Analyse_Var.Results.Cov_fun; %Storing the 'Method' input

% Computing a random normal matrix
RandN = randn(FFTGeometry2D.Ny,FFTGeometry2D.Nx);

% Computing the gaussian field using the spectral representation
Cov = Cov_fun(FFTGeometry2D.X,FFTGeometry2D.Y);
FFTRandN = fft2(RandN);
FFTCov = fft2(fftshift(Cov));
B = real(ifft2(sqrt(FFTCov).*FFTRandN));