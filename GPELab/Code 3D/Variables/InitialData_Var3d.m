%% Computation of the initial data
%% INPUTS:
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var3d.m)
%%          Geometry3D: Structure containing variables concerning the geometry of the problem in 3D (structure) (see Geometry3D_Var3d.m)
%%          Physics3D: Structure containing variables concerning the physics of the problem in 3D (structure) (see Physics3D_Var3d.m)
%% INPUTS(OPTIONAL):
%%          InitialData_Choice: Variable containing the choice between type of computations for the initial data (double, Default: 1)
%%          (Must either be: 1 to compute directly a centered gaussian, 2 to compute directly the Thomas-Fermi approximation or 3 to compute the initial data using the CNSP-CNFG scheme)
%%          X0,Y0,Z0: Coordinates of the center of Gaussians or Thomas-Fermi approximation (vector or double, Default: 0)
%%          gamma_x,gamma_y,gamma_z:Parameters for the centered gaussian (double, Default: 1) (see GaussianInitialData3d.m)
%% OUTPUT:
%%          Phi_0: Initial data computed (cell array)
%% FUNCTIONS USED:
%%          GaussianInitialData3d: To compute the centered gaussian (line 67,84 and 116)
%%          Thomas_Fermi3d: To compute the Thomas-Fermi approximation (line 75 and 103)
%%          CNSP_CNGF3d: To compute initial data using the CNSP-CNFG method (line 110 and 123)

function [Phi_0] = InitialData_Var3d(varargin)
%% Analysis of the inputs
Analyse_Var = inputParser; % Creating the parser
Analyse_Var.addRequired('Method'); % Required input 'Method'
Analyse_Var.addRequired('Geometry3D'); % Required input 'Geometry3D'
Analyse_Var.addRequired('Physics3D'); % Required input 'Physics3D'
Analyse_Var.addOptional('InitialData_Choice',1); % Optional input 'InitialData_Choice' with default value 1
Analyse_Var.addOptional('X0',0); % Optional input 'X0' with default value 0
Analyse_Var.addOptional('Y0',0); % Optional input 'Y0' with default value 0
Analyse_Var.addOptional('Z0',0); % Optional input 'Z0' with default value 0
Analyse_Var.addOptional('gamma_x',1); % Optional input 'gamma_x' with default value 1
Analyse_Var.addOptional('gamma_y',1); % Optional input 'gamma_y' with default value 1
Analyse_Var.addOptional('gamma_z',1); % Optional input 'gamma_z' with default value 1

%% Parsing inputs and storing inputs
% Parsing inputs
Analyse_Var.parse(varargin{:}); % Analysing the inputs
Method = Analyse_Var.Results.Method; %Storing the 'Method' input
Geometry3D = Analyse_Var.Results.Geometry3D; %Storing the 'Geometry2D' input
Physics3D = Analyse_Var.Results.Physics3D; %Storing the 'Physics2D' input
InitialData_Choice = Analyse_Var.Results.InitialData_Choice; %Storing the 'InitialData_Choice' input
X0 = Analyse_Var.Results.X0; %Storing the 'X0' input
% IF X0 is a scalar
if (isscalar(X0))
   X0 = X0*ones(1,Method.Ncomponents); % Storing 'X0' as a vector
end
Y0 = Analyse_Var.Results.Y0; %Storing the 'Y0' input
% IF Y0 is a scalar
if (isscalar(Y0))
   Y0 = Y0*ones(1,Method.Ncomponents); % Storing 'Y0' as a vector
end
Z0 = Analyse_Var.Results.Z0; %Storing the 'Z0' input
% IF Z0 is a scalar
if (isscalar(Z0))
   Z0 = Z0*ones(1,Method.Ncomponents); % Storing 'Z0' as a vector
end
gamma_x = Analyse_Var.Results.gamma_x; %Storing the 'gamma_x' input
gamma_y = Analyse_Var.Results.gamma_y; %Storing the 'gamma_y' input
gamma_z = Analyse_Var.Results.gamma_z; %Storing the 'gamma_y' input


%% Initialization of the initial data variable
Phi_0 = cell(1,Method.Ncomponents); % Initialization of the initial data

%% Computation of a single initial wave function
    % IF one has chosen to compute directly a centered gaussian
    if (InitialData_Choice == 1)
        %FOR each component
        for n = 1:Method.Ncomponents
            Phi_0{n} = GaussianInitialData3d(Geometry3D, Physics3D, gamma_x,gamma_y,gamma_z, X0(n), Y0(n), Z0(n)); % Computing the centered gaussian
            Phi_0{n} = Phi_0{n}/L2_norm3d(Phi_0{n},Geometry3D); %Normalization of the initial function 
        end
    % ELSEIF one has chosen to compute directly the Thomas-Fermi
    % approximation and the nonlinearity is not null
    elseif (InitialData_Choice == 2) && (Physics3D.Beta ~= 0)
        %FOR each component
        for n = 1:Method.Ncomponents
            Phi_0{n} = Thomas_Fermi3d(gamma_x,gamma_y,gamma_z,Physics3D.Beta,Physics3D.Potential_function{n,n}(Geometry3D.X-X0(n),Geometry3D.Y-Y0(n),Geometry3D.Z-Z0(n)) + Physics3D.TimePotential_function{n,n}(0,Geometry3D.X-X0(n),Geometry3D.Y-Y0(n),Geometry3D.Z-Z0(n))); % Computing the Thomas-Fermi approximation
            Phi_0{n} = Phi_0{n}/L2_norm3d(Phi_0{n},Geometry3D); %Normalization of the initial function 
        end
    % ELSEIF one has chosen to compute directly the Thomas-Fermi
    % approximation and the nonlinearity is null
    elseif (InitialData_Choice == 2) && (Physics3D.Beta == 0)
        fprintf('Cannot compute initial data using Thomas_Fermi2d with coefficient Beta = 0. Using the centered gaussian instead.\n');
        %FOR each component
        for n = 1:Method.Ncomponents
            Phi_0{n} = GaussianInitialData3d(Geometry3D, Physics3D, gamma_x,gamma_y, gamma_z, X0(n), Y0(n), Z0(n)); % Computing the centered gaussian
            Phi_0{n} = Phi_0{n}/L2_norm3d(Phi_0{n},Geometry3D); %Normalization of the initial function 
        end
    % ELSEIF one has chosen to compute the initial data using the CNSP-CNFG
    elseif (InitialData_Choice == 3)
        %% Initialization 
        Print.Print = 1; % Not printing informations
        Print.Evo = 5; % Initializing Print.Evo variable
        Print.Draw = 0; % Not drawing functions
        Method.Deltat = 1e-1; % Step time set to 0.1
        Method.Precond = 'ThomasFermi'; % Using the Thomas-Fermi preconditionner
        Outputs = OutputsINI_Var3d(Method);
        Figure = []; % No figure informations
        % IF the nonlinearity is not null
        if (Physics3D.Beta ~= 0)
            fprintf('Computing initial data with the BESP scheme using Thomas_Fermi3d for initial function\n');
            % FOR each component
            for n = 1:Method.Ncomponents
                Phi_0{n} = Thomas_Fermi3d(gamma_x,gamma_y,gamma_z,Physics3D.Beta,Physics3D.Potential_function{n,n}(Geometry3D.X,Geometry3D.Y,Geometry3D.Z)); % Computing the Thomas-Fermi approximation
                % IF the diagonal nonlinearity has been defined
                if (isempty(Physics3D.Nonlinearity_function_Index{n}) == 0)
                    Physics3D.Nonlinearity_function_Index{n} = n; % Restriction to the diagonal nonlinearities
                % ELSE if the diagonal nonlinearity has not been defined
                else
                    Physics3D.Nonlinearity_function_Index{n} = []; % Restriction to the diagonal nonlinearities
                end
                % IF the diagonal non local nonlinearity has been defined
                if (isempty(Physics3D.FFTNonlinearity_function_Index{n}) == 0)
                    Physics3D.FFTNonlinearity_function_Index{n} = n; % Restriction to the diagonal nonlinearities
                % ELSE if the diagonal non local nonlinearity has not been defined
                else
                    Physics3D.FFTNonlinearity_function_Index{n} = []; % Restriction to the diagonal nonlinearities
                end
                % IF the diagonal gradient in the x direction has been defined
                if (isempty(Physics3D.Gradientx_function_Index{n}) == 0)
                    Physics3D.Gradientx_function_Index{n} = n; % Restriction to the diagonal gradients in the x direction
                % ELSE  if the diagonal gradient in the x direction has not been defined
                else
                    Physics3D.Gradientx_function_Index{n} = []; % Restriction to the diagonal gradients in the x direction
                end
                % IF the diagonal gradient in the y direction has been defined
                if (isempty(Physics3D.Gradienty_function_Index{n}) == 0)
                    Physics3D.Gradienty_function_Index{n} = n; % Restriction to the diagonal gradients in the y direction
                % ELSE if the diagonal gradient in the y direction has not been defined
                else
                    Physics3D.Gradienty_function_Index{n} = []; % Restriction to the diagonal gradients in the y direction
                end
                % IF the diagonal gradient in the z direction has been defined
                if (isempty(Physics3D.Gradientz_function_Index{n}) == 0)
                    Physics3D.Gradientz_function_Index{n} = n; % Restriction to the diagonal gradients in the z direction
                % ELSE if the diagonal gradient in the y direction has not been defined
                else
                    Physics3D.Gradientz_function_Index{n} = []; % Restriction to the diagonal gradients in the z direction
                end
                % IF the diagonal potential has been defined
                if (isempty(Physics3D.Potential_function_Index{n}) == 0)
                    Physics3D.Potential_function_Index{n} = n; % Restriction to the diagonal potentials
                % ELSE if the diagonal potential has been defined
                else
                    Physics3D.Potential_function_Index{n} = []; % Restriction to the diagonal potentials
                end
            end
            Phi_0 = BESP_CNGF3d(Phi_0, Method, Geometry3D, Physics3D, Outputs, Print, Figure); % Computing the initial function using the CNSP scheme         
        % ELSEIF the nonlinearity is null
        else
            fprintf('Computing initial data with the BESP scheme using GaussianInitialData3d for initial function\n');
            % FOR each component
            for n = 1:Method.Ncomponents
                Phi_0{n} = GaussianInitialData3d(Geometry3D, Physics3D, gamma_x,gamma_y,gamma_z, X0, Y0, Z0); % Computing the centered gaussian
                % IF the diagonal nonlinearity has been defined
                if (isempty(Physics3D.Nonlinearity_function_Index{n}) == 0)
                    Physics3D.Nonlinearity_function_Index{n} = n; % Restriction to the diagonal nonlinearities
                % ELSE if the diagonal nonlinearity has not been defined
                else
                    Physics3D.Nonlinearity_function_Index{n} = []; % Restriction to the diagonal nonlinearities
                end
                % IF the diagonal non local nonlinearity has been defined
                if (isempty(Physics3D.FFTNonlinearity_function_Index{n}) == 0)
                    Physics3D.FFTNonlinearity_function_Index{n} = n; % Restriction to the diagonal nonlinearities
                % ELSE if the diagonal non local nonlinearity has not been defined
                else
                    Physics3D.FFTNonlinearity_function_Index{n} = []; % Restriction to the diagonal nonlinearities
                end
                % IF the diagonal gradient in the x direction has been defined
                if (isempty(Physics3D.Gradientx_function_Index{n}) == 0)
                    Physics3D.Gradientx_function_Index{n} = n; % Restriction to the diagonal gradients in the x direction
                % ELSE  if the diagonal gradient in the x direction has not been defined
                else
                    Physics3D.Gradientx_function_Index{n} = []; % Restriction to the diagonal gradients in the x direction
                end
                % IF the diagonal gradient in the y direction has been defined
                if (isempty(Physics3D.Gradienty_function_Index{n}) == 0)
                    Physics3D.Gradienty_function_Index{n} = n; % Restriction to the diagonal gradients in the y direction
                % ELSE if the diagonal gradient in the y direction has not been defined
                else
                    Physics3D.Gradienty_function_Index{n} = []; % Restriction to the diagonal gradients in the y direction
                end
                % IF the diagonal gradient in the z direction has been defined
                if (isempty(Physics3D.Gradientz_function_Index{n}) == 0)
                    Physics3D.Gradientz_function_Index{n} = n; % Restriction to the diagonal gradients in the z direction
                % ELSE if the diagonal gradient in the y direction has not been defined
                else
                    Physics3D.Gradientz_function_Index{n} = []; % Restriction to the diagonal gradients in the z direction
                end
                % IF the diagonal potential has been defined
                if (isempty(Physics3D.Potential_function_Index{n}) == 0)
                    Physics3D.Potential_function_Index{n} = n; % Restriction to the diagonal potentials
                % ELSE if the diagonal potential has been defined
                else
                    Physics3D.Potential_function_Index{n} = []; % Restriction to the diagonal potentials
                end
            end
            Phi_0 = BESP_CNGF3d(Phi_0, Method, Geometry3D, Physics3D, Outputs, Print, Figure); % Computing the initial function using the CNSP scheme 
        end
        % FOR each component
        for n = 1:Method.Ncomponents
            Phi_0{n} = Phi_0{n}/L2_norm3d(Phi_0{n},Geometry3D); %Normalization of the initial function   
        end
    end