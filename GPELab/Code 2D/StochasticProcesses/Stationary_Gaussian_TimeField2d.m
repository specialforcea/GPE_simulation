function [B] = Stationary_Gaussian_TimeField2d(Method,Geometry2D,varargin)
%% Analysis of the inputs
Analyse_Var = inputParser; % Creating the parser
Analyse_Var.addOptional('Cov_fun', @(t,X,Y) exp(-X.^2-Y.^2-t^2)); % Required input 'Cov_fun' corresponding to a covariance function

%% Parsing inputs and storing inputs
% Parsing inputs
Analyse_Var.parse(varargin{:}); % Analysing the inputs
Cov_fun = Analyse_Var.Results.Cov_fun; %Storing the 'Method' input

% Computing the time vector
TimeVec = Method.Deltat:Method.Deltat:ceil(Method.Stop_time/Method.Deltat);

% Computing a random normal matrix
RandN = randn(Geometry2D.Ny,Geometry2D.Nx);

% Computing the gaussian field using the spectral representation
Cov = Cov_fun(X,Y);
FFTRandN = fft2(RandN);
FFTCov = fft2(fftshift(Cov));
B = real(ifft2(sqrt(FFTCov).*FFTRandN));