%% Creation of the physiscs structure
%% INPUT:
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var3d.m)
%% INPUTS(OPTIONAL):
%%          Delta: Coefficient in front of the laplacian operator (double)
%%          Beta: Coefficient in front of the nonlinearity (double)
%%          Omega: Coefficient in front of the default gradient function (double)
%% OUTPUT:
%%          Physics3D: Structure containing variables concerning the physics of the problem in 3D (structure)

function [Physics3D] = Physics3D_Var3d(Method, varargin)

%% Physical Inputs
Analyse_Var = inputParser; % Creating the parser
Analyse_Var.addOptional('Delta',1/2,@(x)isscalar(x)); % Optional input 'Delta' with default value 1/2
Analyse_Var.addOptional('Beta',1,@(x)isscalar(x)); % Optional input 'Beta' with default value 0
Analyse_Var.addOptional('Omega',[0,0,0],@(x)isnumeric(x)); % Optional input 'Omega' with default value [0,0,0]
Analyse_Var.parse(varargin{:}); % Analysing the inputs

%% Contructing the physical structure
Physics3D.Delta = Analyse_Var.Results.Delta; % Storing the 'Delta' input
Physics3D.Beta = Analyse_Var.Results.Beta; % Storing the 'Beta' input
% IF Omega is a scalar
if (isscalar(Analyse_Var.Results.Omega))
    Physics3D.Omega = Analyse_Var.Results.Omega*[0,0,1]; % Storing the 'Omega' input as a vector in the z axis
% ELSEIF Omega is a vector of length 3
elseif (isvector(Analyse_Var.Results.Omega)) && (length(Analyse_Var.Results.Omega) == 3)
    Physics3D.Omega = Analyse_Var.Results.Omega; % Storing the 'Omega' input
% ELSE if Omega has not a length of 3
else
    Physics3D.Omega = zeros(1,3); % Initializing the 'Omega' input
    Physics3D.Omega(1:min(length(Analyse_Var.Results.Omega),3)) = Analyse_Var.Results.Omega1:min(length(Analyse_Var.Results.Omega),3); % Storing the 'Omega' input
end

%% Initializing functions
Physics3D.Potential_function = cell(Method.Ncomponents);
Physics3D.Dispersion_function = cell(Method.Ncomponents);
Physics3D.TimeDispersion_function = cell(Method.Ncomponents);
Physics3D.TimePotential_function = cell(Method.Ncomponents);
Physics3D.StochasticPotential_function = cell(Method.Ncomponents);
Physics3D.StochasticDispersion_function = cell(Method.Ncomponents);
Physics3D.IntegratedTimePotential_function = cell(Method.Ncomponents);
Physics3D.IntegratedTimeDispersion_function = cell(Method.Ncomponents);
Physics3D.Nonlinearity_function = cell(Method.Ncomponents);
Physics3D.Nonlinearity_energy_function = cell(Method.Ncomponents);
Physics3D.FFTNonlinearity_function = cell(Method.Ncomponents);
Physics3D.FFTNonlinearity_energy_function = cell(Method.Ncomponents);
Physics3D.Gradientx_function = cell(Method.Ncomponents);
Physics3D.Gradienty_function = cell(Method.Ncomponents);
Physics3D.Gradientz_function = cell(Method.Ncomponents);
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        Physics3D.Potential_function{n,m} = @(X,Y,Z) 0;
        Physics3D.Dispersion_function{n,m} = @(FFTX,FFTY,FFTZ) 0;
        Physics3D.TimeDispersion_function{n,m} = @(t,FFTX,FFTY,FFTZ) 0;
        Physics3D.TimePotential_function{n,m} = @(t,X,Y,Z) 0;
        Physics3D.StochasticPotential_function{n,m} = @(W,X,Y,Z) 0;
        Physics3D.StochasticDispersion_function{n,m} = @(W,FFTX,FFTY,FFTZ) 0;
        Physics3D.StochasticProcess_function = @(t,X,Y,Z) 0;
        Physics3D.Nonlinearity_function{n,m} = @(Phi,X,Y,Z) 0;
        Physics3D.FFTNonlinearity_function{n,m} = @(Phi,X,Y,Z,FFTX,FFTY,FTTZ) 0;
        Physics3D.Gradientx_function{n,m} = @(X,Y,Z) 0;
        Physics3D.Gradienty_function{n,m} = @(X,Y,Z) 0;
        Physics3D.Gradientz_function{n,m} = @(X,Y,Z) 0;
    end
end

%% Initializing index functions
Physics3D.Potential_function_Index = cell(1,Method.Ncomponents);
Physics3D.Dispersion_function_Index = cell(1,Method.Ncomponents);
Physics3D.TimeDispersion_function_Index = cell(1,Method.Ncomponents);
Physics3D.TimePotential_function_Index = cell(1,Method.Ncomponents);
Physics3D.StochasticPotential_function_Index = cell(1,Method.Ncomponents);
Physics3D.StochasticDispersion_function_Index = cell(1,Method.Ncomponents);
Physics3D.Nonlinearity_function_Index = cell(1,Method.Ncomponents);
Physics3D.FFTNonlinearity_function_Index = cell(1,Method.Ncomponents);
Physics3D.Gradientx_function_Index = cell(1,Method.Ncomponents);
Physics3D.Gradienty_function_Index = cell(1,Method.Ncomponents);
Physics3D.Gradientz_function_Index = cell(1,Method.Ncomponents);

%% Initializing gradient to compute index
Physics3D.Gradientx_compute_Index = cell(1,Method.Ncomponents);
Physics3D.Gradienty_compute_Index = cell(1,Method.Ncomponents);
Physics3D.Gradientz_compute_Index = cell(1,Method.Ncomponents);
% FOR each component
for n = 1:Method.Ncomponents
    Physics3D.Gradientx_compute_Index{n} = 0; % Initializing the gradient in the x direction for the component as "not to compute"
    Physics3D.Gradienty_compute_Index{n} = 0; % Initializing the gradient in the y direction for the component as "not to compute"
    Physics3D.Gradientz_compute_Index{n} = 0; % Initializing the gradient in the z direction for the component as "not to compute"
end
