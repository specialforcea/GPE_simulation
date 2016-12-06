%% Creation of the physiscs structure
%% INPUT:
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%% INPUTS(OPTIONAL):
%%          Delta: Coefficient in front of the laplacian operator (double)
%%          Beta: Coefficient in front of the nonlinearity (double)
%%          Omega: Coefficient in front of the default gradient function (double)
%% OUTPUT:
%%          Physics2D: Structure containing variables concerning the physics of the problem in 2D (structure)

function [Physics2D] = Physics2D_Var2d(Method, varargin)

%% Physical Inputs
Analyse_Var = inputParser; % Creating the parser
Analyse_Var.addOptional('Delta',1/2,@(x)isscalar(x)); % Optional input 'Delta' with default value 1/2
Analyse_Var.addOptional('Beta',1,@(x)isscalar(x)); % Optional input 'Beta' with default value 0
Analyse_Var.addOptional('Omega',0,@(x)isscalar(x)); % Optional input 'Omega' with default value 0
Analyse_Var.parse(varargin{:}); % Analysing the inputs

%% Contructing the physical structure
Physics2D.Delta = Analyse_Var.Results.Delta; % Storing the 'Delta' input
Physics2D.Beta = Analyse_Var.Results.Beta; % Storing the 'Beta' input
Physics2D.Omega = Analyse_Var.Results.Omega; % Storing the 'Omega' input

%% Initializing functions
Physics2D.Potential_function = cell(Method.Ncomponents);
Physics2D.Dispersion_function = cell(Method.Ncomponents);
Physics2D.TimeDispersion_function = cell(Method.Ncomponents);
Physics2D.TimePotential_function = cell(Method.Ncomponents);
Physics2D.StochasticPotential_function = cell(Method.Ncomponents);
Physics2D.StochasticDispersion_function = cell(Method.Ncomponents);
Physics2D.IntegratedTimePotential_function = cell(Method.Ncomponents);
Physics2D.IntegratedTimeDispersion_function = cell(Method.Ncomponents);
Physics2D.Nonlinearity_function = cell(Method.Ncomponents);
Physics2D.Nonlinearity_energy_function = cell(Method.Ncomponents);
Physics2D.FFTNonlinearity_function = cell(Method.Ncomponents);
Physics2D.FFTNonlinearity_energy_function = cell(Method.Ncomponents);
Physics2D.Gradientx_function = cell(Method.Ncomponents);
Physics2D.Gradienty_function = cell(Method.Ncomponents);
Physics2D.TimeGradientx_function = cell(Method.Ncomponents);
Physics2D.TimeGradienty_function = cell(Method.Ncomponents);
Physics2D.GradientNLx_function = cell(Method.Ncomponents);
Physics2D.GradientNLy_function = cell(Method.Ncomponents);
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        Physics2D.Potential_function{n,m} = @(X,Y) 0;
        Physics2D.Dispersion_function{n,m} = @(FFTX,FFTY) 0;
        Physics2D.TimeDispersion_function{n,m} = @(t,FFTX,FFTY) 0;
        Physics2D.TimePotential_function{n,m} = @(t,X,Y) 0;
        Physics2D.StochasticPotential_function{n,m} = @(W,X,Y) 0;
        Physics2D.StochasticDispersion_function{n,m} = @(W,FFTX,FFTY) 0;
        Physics2D.StochasticProcess_function = @(t,X,Y) 0;
        Physics2D.Nonlinearity_function{n,m} = @(Phi,X,Y) 0;
        Physics2D.FFTNonlinearity_function{n,m} = @(Phi,X,Y,FFTX,FFTY) 0;
        Physics2D.Gradientx_function{n,m} = @(X,Y) 0;
        Physics2D.Gradienty_function{n,m} = @(X,Y) 0;
        Physics2D.TimeGradientx_function{n,m} = @(t,X,Y) 0;
        Physics2D.TimeGradienty_function{n,m} = @(t,X,Y) 0;
        Physics2D.GradientNLx_function{n,m} = @(Phi,X,Y,FFTX,FFTY) 0;
        Physics2D.GradientNLy_function{n,m} = @(Phi,X,Y,FFTX,FFTY) 0;
    end
end

%% Initializing index functions
Physics2D.Potential_function_Index = cell(1,Method.Ncomponents);
Physics2D.Dispersion_function_Index = cell(1,Method.Ncomponents);
Physics2D.TimeDispersion_function_Index = cell(1,Method.Ncomponents);
Physics2D.TimePotential_function_Index = cell(1,Method.Ncomponents);
Physics2D.StochasticPotential_function_Index = cell(1,Method.Ncomponents);
Physics2D.StochasticDispersion_function_Index = cell(1,Method.Ncomponents);
Physics2D.Nonlinearity_function_Index = cell(1,Method.Ncomponents);
Physics2D.FFTNonlinearity_function_Index = cell(1,Method.Ncomponents);
Physics2D.Gradientx_function_Index = cell(1,Method.Ncomponents);
Physics2D.Gradienty_function_Index = cell(1,Method.Ncomponents);
Physics2D.TimeGradientx_function_Index = cell(1,Method.Ncomponents);
Physics2D.TimeGradienty_function_Index = cell(1,Method.Ncomponents);
Physics2D.GradientNLx_function_Index = cell(1,Method.Ncomponents);
Physics2D.GradientNLy_function_Index = cell(1,Method.Ncomponents);

%% Initializing gradient to compute index
Physics2D.Gradientx_compute_Index = cell(1,Method.Ncomponents);
Physics2D.Gradienty_compute_Index = cell(1,Method.Ncomponents);
Physics2D.TimeGradientx_compute_Index = cell(1,Method.Ncomponents);
Physics2D.TimeGradienty_compute_Index = cell(1,Method.Ncomponents);
Physics2D.GradientNLx_compute_Index = cell(1,Method.Ncomponents);
Physics2D.GradientNLy_compute_Index = cell(1,Method.Ncomponents);
% FOR each component
for n = 1:Method.Ncomponents
    Physics2D.Gradientx_compute_Index{n} = 0; % Initializing the gradient in the x direction for the component as "not to compute"
    Physics2D.Gradienty_compute_Index{n} = 0; % Initializing the gradient in the y direction for the component as "not to compute"
    Physics2D.TimeGradientx_compute_Index{n} = 0; % Initializing the gradient in the x direction for the component as "not to compute"
    Physics2D.TimeGradienty_compute_Index{n} = 0; % Initializing the gradient in the y direction for the component as "not to compute"
    Physics2D.GradientNLx_compute_Index{n} = 0; % Initializing the gradient in the x direction for the component as "not to compute"
    Physics2D.GradientNLy_compute_Index{n} = 0; % Initializing the gradient in the y direction for the component as "not to compute"
end
