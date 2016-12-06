function [B] = Stationary_Gaussian_Process2d(Method,varargin)
%% Analysis of the inputs
Analyse_Var = inputParser; % Creating the parser
Analyse_Var.addOptional('Cov_fun', @(t) exp(-t.^2)); % Required input 'Cov_fun' corresponding to a covariance function

%% Parsing inputs and storing inputs
% Parsing inputs
Analyse_Var.parse(varargin{:}); % Analysing the inputs
Cov_fun = Analyse_Var.Results.Cov_fun; %Storing the 'Method' input

% Computing the time vector
TimeVec = Method.Deltat:Method.Deltat:Method.Stop_time;

% Computing a random normal vector
RandN = randn(1,length(TimeVec));

% Computing the gaussian process using the spectral representation
Cov = Cov_fun(TimeVec);
FFTRandN = fft(RandN);
FFTCov = fft(fftshift(Cov));
BVect = real(ifft(sqrt(FFTCov).*FFTRandN));
B = @(t) BVect(t/Method.Deltat);