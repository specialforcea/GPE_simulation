%% Launching simulation and applying the continuation method
%% INPUTS:
%%          Phi_0: Initial wave functions (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var3d.m)
%%          Geometry3D: Structure containing variables concerning the geometry of the problem in 3D (structure) (see Geometry3D_Var3d.m)
%%          Physics3D: Structure containing variables concerning the physics of the problem in 3D (structure) (see Physics3D_Var3d.m)
%%          Outputs: Different outputs computed during the computation of the ground states (structure) (see OutputsINI_Var3D.m and FFTOutputs_Var3d.m)
%%          Print: Structure containing variables concerning the printing and drawing of informations during the computations (structure) (see Print_Var3d.m)
%%          Figure: Structure containing variables concerning the figures (structure)
%% INPUTS(OPTIONAL):
%%          Continuation: Continuation structure used to apply a continuation method (structure) (see Continuation_Var3d.m)
%% OUTPUTS:
%%          Phi: Wave function at the end of the simulation (cell array)
%%          Outputs: Different outputs computed during the simulation (structure) (see OutputsINI_Var3d.m and FFTOutputs_Var3d.m)
%% FUNCTIONS USED:
%%          Physics3D_Var3d: To replace the Delta, Beta or Omega coefficients for the continuation method (line 89,92 and 95)
%%          Potential_Var3d: To replace the G matrix of the potential for the continuation method (line 98)
%%          Nonlinearity_Var3d: To replace the G matrix of the nonlinearity for the continuation method (line 101)
%%          Gradientx_Var3d: To replace the G matrix of the gradientx for the continuation method (line 104)
%%          Gradienty_Var3d: To replace the G matrix of the gradienty for the continuation method (line 107)
%%          Gradientz_Var3d: To replace the G matrix of the gradienty for the continuation method (line 110)
%%          BESP_CNGF3d: To compute a ground state using the BESP scheme (line 117)
%%          CNSP_CNGF3d: To compute a ground state using the CNSP scheme (line 119)
%%          BEFD_CNGF3d: To compute a ground state using the BEFD scheme (line 121)
%%          CNFD_CNGF3d: To compute a ground state using the CNFD scheme (line 123)
%%          BESP_DS3d: To compute a dynamic solution using the BESP scheme (line 128)
%%          CNSP_DS3d: To compute a dynamic solution using the CNSP scheme (line 130)
%%          BEFD_DS3d: To compute a dynamic solution using the BEFD scheme (line 132)
%%          CNFD_DS3d: To compute a dynamic solution using the CNFD scheme (line 134)

function [Phi, Outputs] = GPELab3d(varargin)
%% Analysis of the inputs
Analyse_Var = inputParser; % Creating the parser
Analyse_Var.addRequired('Phi_0',@(x)iscell(x)); % Required input 'Phi_0' which must be a matrix
Analyse_Var.addRequired('Method',@(x)isstruct(x)); % Required input 'Method' which must be a structure
Analyse_Var.addRequired('Geometry3D',@(x)isstruct(x)); % Required input 'Geometry3D' which must be a structure
Analyse_Var.addRequired('Physics3D',@(x)isstruct(x)); % Required input 'Physics3D' which must be a structure
Analyse_Var.addRequired('Outputs',@(x)isstruct(x)); % Required input 'Outputs' which must be a structure
Analyse_Var.addOptional('Continuation',[],@(x)isstruct(x)+isempty(x)); % Optional input 'Continuation' which must be a structure
Analyse_Var.addOptional('Print', Print_Var3d(1,10,1),@(x)isstruct(x)); % Optional input 'Print' which must be a structure
Analyse_Var.addOptional('Figure', Figure_Var3d,@(x)isstruct(x)); % Optional input 'Figure' which must be a structure

%% Parsing inputs and creating the Method structure
% Parsing inputs
Analyse_Var.parse(varargin{:}); % Analysing the inputs
% Contructing the Method structure
Phi_0 = Analyse_Var.Results.Phi_0; %Storing the 'Phi_0' input
Method = Analyse_Var.Results.Method; %Storing the 'Method' input
Geometry3D = Analyse_Var.Results.Geometry3D; %Storing the 'Geometry3D' input
Physics3D = Analyse_Var.Results.Physics3D; %Storing the 'Physics3D' input
Outputs = Analyse_Var.Results.Outputs; %Storing the 'Outputs' input
Continuation = Analyse_Var.Results.Continuation; %Storing the 'Continuation' input
Print = Analyse_Var.Results.Print; %Storing the 'Print' input
Figure = Analyse_Var.Results.Figure; %Storing the 'Figure' input

% IF Evo_outputs is bigger than Evo
if (Outputs.Evo_outputs > Print.Evo)
    Outputs.Evo_outputs = Print.Evo; % Setting the two variables as the same
end
% IF the continuation structure is provided
if (isempty(Continuation) == 0)
    % FOR each coefficient
    for n = 1:length(Continuation.Coefficient)
        % IF it is the first coefficient
        if (n == 1)
            Continuation_iter = length(Continuation.Coefficient{n}); % Storing the number of iteration for the continuation as the length of the coefficient evolution
        % ELSE if it is not the first coefficient
        else
            Continuation_iter = min(Continuation_iter,length(Continuation.Coefficient{n})); % Storing the number of iteration for the continuation as minimum between the length of the coefficient evolution or the previous number
        end
    end
% Storing temporary functions
TMP_Potential_function = Physics3D.Potential_function; % Storing temporary potential functions
TMP_Nonlinearity_function = Physics3D.Nonlinearity_function; % Storing temporary nonlinearity functions
TMP_Nonlinearity_energy_function = Physics3D.Nonlinearity_energy_function; % Storing temporary nonlinearity functions
TMP_FFTNonlinearity_function = Physics3D.FFTNonlinearity_function; % Storing temporary nonlinearity functions
TMP_FFTNonlinearity_energy_function = Physics3D.FFTNonlinearity_energy_function; % Storing temporary non-local nonlinearity functions
TMP_Gradientx_function = Physics3D.Gradientx_function; % Storing temporary gradientx functions
TMP_Gradienty_function = Physics3D.Gradienty_function; % Storing temporary gradienty functions
TMP_Gradientz_function = Physics3D.Gradientz_function; % Storing temporary gradienty functions
% ELSE if the continuation structure is not provided
else
    Continuation_iter = 1; % Storing the number of iteration for the continuation as 1
end

%% Continuation
for n = 1:Continuation_iter
    % IF the continuation structure is provided
    if (isempty(Continuation) == 0)
       % FOR each coefficient that needs to be modified during the
       % continuation
       for m = 1:Continuation.Ncontinuation
           % IF the coefficient is Delta
           if (strcmp(Continuation.Coefficient_name{m},'Delta'))
               Physics3D = Physics3D_Var3d(Method, Continuation.Coefficient{m}{n},Physics3D.Beta,Physics3D.Omega); % Replace the coefficient Delta
           % ELSEIF the coefficient is Beta
           elseif (strcmp(Continuation.Coefficient_name{m},'Beta'))
               Physics3D = Physics3D_Var3d(Method,Physics3D.Delta,Continuation.Coefficient{m}{n},Physics3D.Omega); % Replace the coefficient Beta
           % ELSEIF the coefficient is Omega
           elseif (strcmp(Continuation.Coefficient_name{m},'Omega'))
               Physics3D = Physics3D_Var3d(Method,Physics3D.Delta,Physics3D.Beta,Continuation.Coefficient{m}{n}); % Replace the coefficient Omega
           % ELSEIF the coefficient is the G matrix of the potential
           elseif (strcmp(Continuation.Coefficient_name{m},'GPotential'))
               Physics3D = Potential_Var3d(Method, Physics3D, TMP_Potential_function, Continuation.Coefficient{m}{n}); % Replace the G matrix of the potential
           % ELSEIF the coefficient is the G matrix of the nonlinearity
           elseif (strcmp(Continuation.Coefficient_name{m},'GNonlinearity'))
               Physics3D = Nonlinearity_Var3d(Method, Physics3D, TMP_Nonlinearity_function, Continuation.Coefficient{m}{n}, TMP_Nonlinearity_energy_function); % Replace the G matrix of the nonlinearity
           % ELSEIF the coefficient is the G matrix of the nonlinearity
           elseif (strcmp(Continuation.Coefficient_name{m},'GFFTNonlinearity'))
               Physics3D = FFTNonlinearity_Var3d(Method, Physics3D, TMP_FFTNonlinearity_function, Continuation.Coefficient{m}{n}, TMP_FFTNonlinearity_energy_function); % Replace the G matrix of the non-local nonlinearity
           % ELSEIF the coefficient is the G matrix of the gradientx
           elseif (strcmp(Continuation.Coefficient_name{m},'GGradientx'))
               Physics3D = Gradientx_Var3d(Method, Physics3D, TMP_Gradientx_function, Continuation.Coefficient{m}{n}); % Replace the G matrix of the gradientx
           % ELSEIF the coefficient is the G matrix of the gradienty
           elseif (strcmp(Continuation.Coefficient_name{m},'GGradienty'))
               Physics3D = Gradienty_Var3d(Method, Physics3D, TMP_Gradienty_function, Continuation.Coefficient{m}{n}); % Replace the G matrix of the gradienty
           % ELSEIF the coefficient is the G matrix of the gradientz
           elseif (strcmp(Continuation.Coefficient_name{m},'GGradientz'))
               Physics3D = Gradientz_Var3d(Method, Physics3D, TMP_Gradientz_function, Continuation.Coefficient{m}{n}); % Replace the G matrix of the gradientz
           end
       end
    end
    %% Ground state computation
    if (strcmp(Method.Computation,'Ground'))
        if (strcmp(Method.Type,'BESP'))
        [Phi,Outputs] = BESP_CNGF3d(Phi_0, Method, Geometry3D, Physics3D, Outputs, Print, Figure);
        elseif (strcmp(Method.Type,'CNSP'))
        [Phi,Outputs] = CNSP_CNGF3d(Phi_0, Method, Geometry3D, Physics3D, Outputs, Print, Figure);
        elseif (strcmp(Method.Type,'BEFD'))
        [Phi, Outputs] = BEFD_CNGF3d(Phi_0, Method, Geometry3D, Physics3D, Outputs, Print, Figure);
        elseif (strcmp(Method.Type,'CNFD'))
        [Phi, Outputs] = CNFD_CNGF3d(Phi_0, Method, Geometry3D, Physics3D, Outputs, Print, Figure);
        end
    %% Dynamic solution computation
    elseif (strcmp(Method.Computation,'Dynamic'))
        if (strcmp(Method.Type,'BESP'))
        [Phi,Outputs] = BESP_DS3d(Phi_0, Method, Geometry3D, Physics3D, Outputs, Print, Figure);
        elseif (strcmp(Method.Type,'CNSP'))
        [Phi,Outputs] = CNSP_DS3d(Phi_0, Method, Geometry3D, Physics3D, Outputs, Print, Figure);
        elseif (strcmp(Method.Type,'BEFD'))
        [Phi, Outputs] = BEFD_DS3d(Phi_0, Method, Geometry3D, Physics3D, Outputs, Print, Figure);
        elseif (strcmp(Method.Type,'CNFD'))
        [Phi, Outputs] = CNFD_DS3d(Phi_0, Method, Geometry3D, Physics3D, Outputs, Print, Figure);
        elseif (strcmp(Method.Type,'Relaxation'))
        [Phi, Outputs] = RSP_DS3d(Phi_0, Method, Geometry3D, Physics3D, Outputs, Print, Figure);
        elseif (strcmp(Method.Type,'Splitting'))
        [Phi, Outputs] = SPL_DS3d(Phi_0, Method, Geometry3D, Physics3D, Outputs, Print, Figure);
        end
    end

end