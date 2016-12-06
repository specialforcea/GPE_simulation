%% Creation of the physiscs structure
%% INPUT:
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var1d.m)
%% INPUTS(OPTIONAL):
%%          Delta: Coefficient in front of the laplacian operator (double)
%%          Beta: Coefficient in front of the nonlinearity (double)
%% OUTPUT:
%%          Physics1D: Structure containing variables concerning the physics of the problem in 1D (structure)

function [Physics1D] = Physics1D_Var1d(Method, varargin)

%% Physical Inputs
Analyse_Var = inputParser; % Creating the parser
Analyse_Var.addOptional('Delta',1/2,@(x)isscalar(x)); % Optional input 'Delta' with default value 1/2
Analyse_Var.addOptional('Beta',1,@(x)isscalar(x)); % Optional input 'Beta' with default value 0
Analyse_Var.parse(varargin{:}); % Analysing the inputs

%% Contructing the physical structure
Physics1D.Delta = Analyse_Var.Results.Delta; % Storing the 'Delta' input
Physics1D.Beta = Analyse_Var.Results.Beta; % Storing the 'Beta' input

%% Initializing functions
Physics1D.Potential_function = cell(Method.Ncomponents);
Physics1D.Disperion_function = cell(Method.Ncomponents);
Physics1D.TimeDisperion_function = cell(Method.Ncomponents);
Physics1D.TimePotential_function = cell(Method.Ncomponents);
Physics1D.StochasticPotential_function = cell(Method.Ncomponents);
Physics1D.IntegratedTimePotential_function = cell(Method.Ncomponents);
Physics1D.IntegratedTimeDispersion_function = cell(Method.Ncomponents);
Physics1D.Nonlinearity_function = cell(Method.Ncomponents);
Physics1D.Nonlinearity_energy_function = cell(Method.Ncomponents);
Physics1D.FFTNonlinearity_function = cell(Method.Ncomponents);
Physics1D.FFTNonlinearity_energy_function = cell(Method.Ncomponents);
Physics1D.Gradientx_function = cell(Method.Ncomponents);
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        Physics1D.Potential_function{n,m} = @(X) 0;
        Physics1D.Dispersion_function{n,m} = @(FFTX) 0;
        Physics1D.TimeDispersion_function{n,m} = @(t,FFTX) 0;
        Physics1D.TimePotential_function{n,m} = @(t,X) 0;
        Physics1D.StochasticPotential_function{n,m} = @(W,X) 0;
        Physics1D.StochasticDispersion_function{n,m} = @(W,FFTX) 0;
        Physics1D.StochasticProcess_function = @(t,X) 0;
        Physics1D.Nonlinearity_function{n,m} = @(Phi,X) 0;
        Physics1D.FFTNonlinearity_function{n,m} = @(Phi,X,FFTX) 0;
        Physics1D.Gradientx_function{n,m} = @(X) 0;
    end
end

%% Initializing index functions
Physics1D.Potential_function_Index = cell(1,Method.Ncomponents);
Physics1D.Dispersion_function_Index = cell(1,Method.Ncomponents);
Physics1D.TimeDispersion_function_Index = cell(1,Method.Ncomponents);
Physics1D.TimePotential_function_Index = cell(1,Method.Ncomponents);
Physics1D.StochasticPotential_function_Index = cell(1,Method.Ncomponents);
Physics1D.StochasticDispersion_function_Index = cell(1,Method.Ncomponents);
Physics1D.Nonlinearity_function_Index = cell(1,Method.Ncomponents);
Physics1D.FFTNonlinearity_function_Index = cell(1,Method.Ncomponents);
Physics1D.Gradientx_function_Index = cell(1,Method.Ncomponents);

%% Initializing gradient to compute index
Physics1D.Gradientx_compute_Index = cell(1,Method.Ncomponents);
% FOR each component
for n = 1:Method.Ncomponents
    Physics1D.Gradientx_compute_Index{n} = 0; % Initializing the gradient in the x direction for the component as "not to compute"
end
