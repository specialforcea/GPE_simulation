function [G] = Gaussian_Process2d(Method,varargin)
%% Analysis of the inputs
Analyse_Var = inputParser; % Creating the parser
Analyse_Var.addOptional('Mean_fun', @(t) 0); % Required input 'Mean_fun' corresponding to a covariance function
Analyse_Var.addOptional('Covariance_fun', @(t,s) min(t,s)); % Required input 'Covariance_fun' corresponding to a covariance function

%% Parsing inputs and storing inputs
% Parsing inputs
Analyse_Var.parse(varargin{:}); % Analysing the inputs
Mean_fun = Analyse_Var.Results.Mean_fun; %Storing the 'Mean_fun' input
Covariance_fun = Analyse_Var.Results.Covariance_fun; %Storing the 'Covariance_fun' input

% Computing the time vector
TimeVec = Method.Deltat:Method.Deltat:ceil(Method.Stop_time/Method.Deltat);

% Computing the time matrices
[TimeMat_t,TimeMat_s] = meshgrid(TimeVec,TimeVec);

% Computing a random normal vector
RandN = randn(1,length(TimeVec));

% Computing the mean vector and the covariance matrix
Mean = Mean_fun(TimeVec);
Covariance = Covariance_fun(TimeMat_t,TimeMat_s);

% Computing the square root of the covariance matrix
SqrtCovariance = chol(Covariance,'lower');

% Computing the gaussian process
G = Mean + SqrtCovariance*RandN;