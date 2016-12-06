%% Computation of the initial data
%% INPUTS:
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var1d.m)
%%          Geometry1D: Structure containing variables concerning the geometry of the problem in 1D (structure) (see Geometry1D_Var1d.m)
%%          Physics1D: Structure containing variables concerning the physics of the problem in 1D (structure) (see Physics1D_Var1d.m)
%% INPUTS(OPTIONAL):
%%          InitialData_Choice: Variable containing the choice between type of computations for the initial data (double, Default: 1)
%%          (Must either be: 1 to compute directly a centered gaussian, 2 to compute directly the Thomas-Fermi approximation or 3 to compute the initial data using the CNSP-CNFG scheme)
%%          X0: Coordinates of the center of Gaussians or Thomas-Fermi approximation (vector or double, Default: 0)
%%          gamma_x:Parameters for the centered gaussian (double, Default: 1) (see GaussianInitialData1d.m)
%% OUTPUT:
%%          Phi_0: Initial data computed (cell array)
%% FUNCTIONS USED:
%%          GaussianInitialData1d: To compute the centered gaussian (line 49,66 and 96)
%%          Thomas_Fermi1d: To compute the Thomas-Fermi approximation (line 57 and 85)
%%          CNSP_CNGF1d: To compute initial data using the CNSP-CNFG method (line 90 and 101)

function [Phi_0] = InitialData_Var1d(varargin)
%% Analysis of the inputs
Analyse_Var = inputParser; % Creating the parser
Analyse_Var.addRequired('Method'); % Required input 'Method'
Analyse_Var.addRequired('Geometry1D'); % Required input 'Geometry1D'
Analyse_Var.addRequired('Physics1D'); % Required input 'Physics1D'
Analyse_Var.addOptional('InitialData_Choice',1); % Optional input 'InitialData_Choice' with default value 1
Analyse_Var.addOptional('X0', 0); % Optional input 'X0' with default value 0
Analyse_Var.addOptional('gamma_x', 1); % Optional input 'gamma_x' with default value 1

%% Parsing inputs and storing inputs
% Parsing inputs
Analyse_Var.parse(varargin{:}); % Analysing the inputs
Method = Analyse_Var.Results.Method; %Storing the 'Method' input
Geometry1D = Analyse_Var.Results.Geometry1D; %Storing the 'Geometry1D' input
Physics1D = Analyse_Var.Results.Physics1D; %Storing the 'Physics1D' input
InitialData_Choice = Analyse_Var.Results.InitialData_Choice; %Storing the 'InitialData_Choice' input
X0 = Analyse_Var.Results.X0; %Storing the 'X0' input
% IF X0 is a scalar
if (isscalar(X0))
   X0 = X0*ones(1,Method.Ncomponents); % Storing 'X0' as a vector
end
gamma_x = Analyse_Var.Results.gamma_x; %Storing the 'gamma_x' input

%% Initialization of the initial data variable
Phi_0 = cell(1,Method.Ncomponents); % Initialization of the initial data

%% Computation of a single initial wave function
    % IF one has chosen to compute directly a centered gaussian
    if (InitialData_Choice == 1)
        %FOR each component
        for n = 1:Method.Ncomponents
            Phi_0{n} = GaussianInitialData1d(Geometry1D, gamma_x, X0(n)); % Computing the centered gaussian
            Phi_0{n} = Phi_0{n}/L2_norm1d(Phi_0{n},Geometry1D); %Normalization of the initial function 
        end
    % ELSEIF one has chosen to compute directly the Thomas-Fermi
    % approximation and the nonlinearity is not null
    elseif (InitialData_Choice == 2) && (Physics1D.Beta ~= 0)
        %FOR each component
        for n = 1:Method.Ncomponents
            Phi_0{n} = Thomas_Fermi1d(gamma_x,Physics1D.Beta,Physics1D.Potential_function{n,n}(Geometry1D.X-X0(n)) + Physics1D.TimePotential_function{n,n}(0,Geometry1D.X-X0(n))); % Computing the Thomas-Fermi approximation
            Phi_0{n} = Phi_0{n}/L2_norm1d(Phi_0{n},Geometry1D); %Normalization of the initial function 
        end
    % ELSEIF one has chosen to compute directly the Thomas-Fermi
    % approximation and the nonlinearity is null
    elseif (InitialData_Choice == 2) && (Physics1D.Beta == 0)
        fprintf('Cannot compute initial data using Thomas_Fermi1d with coefficient Beta = 0. Using the centered gaussian instead.\n');
        %FOR each component
        for n = 1:Method.Ncomponents
            Phi_0{n} = GaussianInitialData1d(Geometry1D, gamma_x, X0(n)); % Computing the centered gaussian
            Phi_0{n} = Phi_0{n}/L2_norm1d(Phi_0{n},Geometry1D); %Normalization of the initial function 
        end
    % ELSEIF one has chosen to compute the initial data using the CNSP-CNFG
    elseif (InitialData_Choice == 3)
        %% Initialization 
        Print.Print = 0; % Not printing informations
        Print.Type = 2; % Using the CNSP scheme
        Print.Evo = 5; % Initializing Print.Evo variable
        Print.Draw = 0; % Not drawing functions
        Method.Deltat = 1e-1; % Step time set to 0.1
        Method.Precond = 'ThomasFermi'; % Using the Thomas-Fermi preconditionner
        Outputs = OutputsINI_Var1d(Method);
        Figure = []; % No figure informations
        % IF the nonlinearity is not null
        if (Physics1D.Beta ~= 0)
            fprintf('Computing initial data with the BESP scheme using Thomas_Fermi1d for initial function\n');
            % FOR each component
            for n = 1:Method.Ncomponents
                Phi_0{n} = Thomas_Fermi1d(gamma_x,Physics1D.Beta,Physics1D.Potential_function{n,n}(Geometry1D.X)); % Computing the Thomas-Fermi approximation
                % IF the diagonal nonlinearity has been defined
                if (isempty(Physics1D.Nonlinearity_function_Index{n}) == 0)
                    Physics1D.Nonlinearity_function_Index{n} = n; % Restriction to the diagonal nonlinearities
                % ELSE if the diagonal nonlinearity has not been defined
                else
                    Physics1D.Nonlinearity_function_Index{n} = []; % Restriction to the diagonal nonlinearities
                end
                % IF the diagonal non local nonlinearity has been defined
                if (isempty(Physics1D.FFTNonlinearity_function_Index{n}) == 0)
                    Physics1D.FFTNonlinearity_function_Index{n} = n; % Restriction to the diagonal nonlinearities
                % ELSE if the diagonal non local nonlinearity has not been defined
                else
                    Physics1D.FFTNonlinearity_function_Index{n} = []; % Restriction to the diagonal nonlinearities
                end
                % IF the diagonal gradient in the x direction has been defined
                if (isempty(Physics1D.Gradientx_function_Index{n}) == 0)
                    Physics1D.Gradientx_function_Index{n} = n; % Restriction to the diagonal gradients in the x direction
                % ELSE  if the diagonal gradient in the x direction has not been defined
                else
                    Physics1D.Gradientx_function_Index{n} = []; % Restriction to the diagonal gradients in the x direction
                end
                % IF the diagonal potential has been defined
                if (isempty(Physics1D.Potential_function_Index{n}) == 0)
                    Physics1D.Potential_function_Index{n} = n; % Restriction to the diagonal potentials
                % ELSE if the diagonal potential has been defined
                else
                    Physics1D.Potential_function_Index{n} = []; % Restriction to the diagonal potentials
                end
            end
            Phi_0 = BESP_CNGF1d(Phi_0, Method, Geometry1D, Physics1D, Outputs, Print, Figure); % Computing the initial function using the CNSP scheme         
        % ELSEIF the nonlinearity is null
        else
            fprintf('Computing initial data with the BESP scheme using GaussianInitialData2d for initial function\n');
            % FOR each component
            for n = 1:Method.Ncomponents
                Phi_0{n} = GaussianInitialData1d(Geometry1D, gamma_x, X0); % Computing the centered gaussian
                % IF the diagonal nonlinearity has been defined
                if (isempty(Physics1D.Nonlinearity_function_Index{n}) == 0)
                    Physics1D.Nonlinearity_function_Index{n} = n; % Restriction to the diagonal nonlinearities
                % ELSE if the diagonal nonlinearity has not been defined
                else
                    Physics1D.Nonlinearity_function_Index{n} = []; % Restriction to the diagonal nonlinearities
                end
                % IF the diagonal non local nonlinearity has been defined
                if (isempty(Physics1D.FFTNonlinearity_function_Index{n}) == 0)
                    Physics1D.FFTNonlinearity_function_Index{n} = n; % Restriction to the diagonal nonlinearities
                % ELSE if the diagonal non local nonlinearity has not been defined
                else
                    Physics1D.FFTNonlinearity_function_Index{n} = []; % Restriction to the diagonal nonlinearities
                end
                % IF the diagonal gradient in the x direction has been defined
                if (isempty(Physics1D.Gradientx_function_Index{n}) == 0)
                    Physics1D.Gradientx_function_Index{n} = n; % Restriction to the diagonal gradients in the x direction
                % ELSE  if the diagonal gradient in the x direction has not been defined
                else
                    Physics1D.Gradientx_function_Index{n} = []; % Restriction to the diagonal gradients in the x direction
                end
                % IF the diagonal potential has been defined
                if (isempty(Physics1D.Potential_function_Index{n}) == 0)
                    Physics1D.Potential_function_Index{n} = n; % Restriction to the diagonal potentials
                % ELSE if the diagonal potential has been defined
                else
                    Physics1D.Potential_function_Index{n} = []; % Restriction to the diagonal potentials
                end
            end
            Phi_0 = BESP_CNGF1d(Phi_0, Method, Geometry1D, Physics1D, Outputs, Print, Figure); % Computing the initial function using the CNSP scheme 
        end
        % FOR each component
        for n = 1:Method.Ncomponents
            Phi_0{n} = Phi_0{n}/L2_norm1d(Phi_0{n},Geometry1D); %Normalization of the initial function   
        end
    end