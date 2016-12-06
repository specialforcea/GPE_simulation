function [B,W] = Brownian_Process2d(Method,varargin)
%% Analysis of the inputs
Analyse_Var = inputParser; % Creating the parser
Analyse_Var.addOptional('H', 1/2); % Required input 'H' corresponding to a covariance function

%% Parsing inputs and storing inputs
% Parsing inputs
Analyse_Var.parse(varargin{:}); % Analysing the inputs
H = Analyse_Var.Results.H; %Storing the 'H' input

%% Initializing
% Computing the time vector
LengthTime = ceil(Method.Stop_time/Method.Deltat);
TimeVec = [0:LengthTime+1];

% Computing a complex random normal vector
RandN = randn(1,2*LengthTime) + 1i*randn(1,2*LengthTime);

%% Computing the brownian motion using the spectral representation
% Computing the spectrum of the covariance
Covariance = abs(TimeVec + 1).^(2*H)+abs(TimeVec - 1).^(2*H)-2*abs(TimeVec).^(2*H);
SpectCovariance = [Covariance,Covariance(LengthTime-1:-1:2)];
FFTCov = fft(SpectCovariance);
% Computing the increments of the brownian motion
WVect = ifft(sqrt(FFTCov).*RandN)*(Method.Deltat)^(H-1/2);
WVect = real(WVect(1:LengthTime));
% Computing the brownian motion
BVect = [0,cumsum(WVect)];
WVect = [0,WVect];
B = @(t) BVect(round(t/Method.Deltat)+1);
W = @(t) WVect(round(t/Method.Deltat)+1);