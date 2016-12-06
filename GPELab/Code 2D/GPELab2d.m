%% Launching simulation and applying the continuation method
%% INPUTS:
%%          Phi_0: Initial wave functions (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          Geometry2D: Structure containing variables concerning the geometry of the problem in 2D (structure) (see Geometry2D_Var2d.m)
%%          Physics2D: Structure containing variables concerning the physics of the problem in 2D (structure) (see Physics2D_Var2d.m)
%%          Outputs: Different outputs computed during the computation of the ground states (structure) (see OutputsINI_Var2D.m and FFTOutputs_Var2d.m)
%% INPUTS(OPTIONAL):
%%          Continuation: Continuation structure used to apply a continuation method (structure) (see Continuation_Var2d.m)
%%          Print: Structure containing variables concerning the printing and drawing of informations during the computations (structure) (see Print_Var2d.m)
%%          Figure: Structure containing variables concerning the figures (structure)
%% OUTPUTS:
%%          Phi: Wave function at the end of the simulation (cell array)
%%          Outputs: Different outputs computed during the simulation (structure) (see OutputsINI_Var2D.m and FFTOutputs_Var2d.m)
%% FUNCTIONS USED:
%%          Physics2D_Var2d: To replace the Delta, Beta or Omega coefficients for the continuation method (line 87,90 and 93)
%%          Potential_Var2d: To replace the G matrix of the potential for the continuation method (line 96)
%%          Nonlinearity_Var2d: To replace the G matrix of the nonlinearity for the continuation method (line 99)
%%          Gradientx_Var2d: To replace the G matrix of the gradientx for the continuation method (line 102)
%%          Gradienty_Var2d: To replace the G matrix of the gradienty for the continuation method (line 105)
%%          BESP_CNGF2d: To compute a ground state using the BESP scheme (line 112)
%%          CNSP_CNGF2d: To compute a ground state using the CNSP scheme (line 114)
%%          BEFD_CNGF2d: To compute a ground state using the BEFD scheme (line 116)
%%          CNFD_CNGF2d: To compute a ground state using the CNFD scheme (line 118)
%%          BESP_DS2d: To compute a dynamic solution using the BESP scheme (line 123)
%%          CNSP_DS2d: To compute a dynamic solution using the CNSP scheme (line 125)
%%          BEFD_DS2d: To compute a dynamic solution using the BEFD scheme (line 127)
%%          CNFD_DS2d: To compute a dynamic solution using the CNFD scheme (line 129)

function [Phi, Outputs] = GPELab2d(varargin)
%% Analysis of the inputs
Analyse_Var = inputParser; % Creating the parser
Analyse_Var.addRequired('Phi_0',@(x)iscell(x)); % Required input 'Phi_0' which must be a matrix
Analyse_Var.addRequired('Method',@(x)isstruct(x)); % Required input 'Method' which must be a structure
Analyse_Var.addRequired('Geometry2D',@(x)isstruct(x)); % Required input 'Geometry2D' which must be a structure
Analyse_Var.addRequired('Physics2D',@(x)isstruct(x)); % Required input 'Physics2D' which must be a structure
Analyse_Var.addRequired('Outputs',@(x)isstruct(x)); % Required input 'Outputs' which must be a structure
Analyse_Var.addOptional('Continuation',[],@(x)isstruct(x)+isempty(x)); % Optional input 'Continuation' which must be a structure
Analyse_Var.addOptional('Print', Print_Var2d,@(x)isstruct(x)); % Optional input 'Print' which must be a structure
Analyse_Var.addOptional('Figure', Figure_Var2d,@(x)isstruct(x)); % Optional input 'Figure' which must be a structure


%% Parsing inputs and creating the Method structure
% Parsing inputs
Analyse_Var.parse(varargin{:}); % Analysing the inputs
% Contructing the Method structure
Phi_0 = Analyse_Var.Results.Phi_0; %Storing the 'Phi_0' input
Method = Analyse_Var.Results.Method; %Storing the 'Method' input
Geometry2D = Analyse_Var.Results.Geometry2D; %Storing the 'Geometry2D' input
Physics2D = Analyse_Var.Results.Physics2D; %Storing the 'Physics2D' input
Outputs = Analyse_Var.Results.Outputs; %Storing the 'Outputs' input
Continuation = Analyse_Var.Results.Continuation; %Storing the 'Continuation' input
Print = Analyse_Var.Results.Print; %Storing the 'Print' input
Figure = Analyse_Var.Results.Figure; %Storing the 'Figure' input

% IF Evo_outputs is bigger than Evo
if (Outputs.Evo_outputs>Print.Evo)
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
TMP_Potential_function = Physics2D.Potential_function; % Storing temporary potential functions
TMP_Nonlinearity_function = Physics2D.Nonlinearity_function; % Storing temporary nonlinearity functions
TMP_Nonlinearity_energy_function = Physics2D.Nonlinearity_energy_function; % Storing temporary nonlinearity functions
TMP_FFTNonlinearity_function = Physics2D.FFTNonlinearity_function; % Storing temporary non-local nonlinearity functions
TMP_FFTNonlinearity_energy_function = Physics2D.FFTNonlinearity_energy_function; % Storing temporary non-local nonlinearity functions
TMP_Gradientx_function = Physics2D.Gradientx_function; % Storing temporary gradientx functions
TMP_Gradienty_function = Physics2D.Gradienty_function; % Storing temporary gradienty functions
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
               Physics2D = Physics2D_Var2d(Method, Continuation.Coefficient{m}{n},Physics2D.Beta,Physics2D.Omega); % Replace the coefficient Delta
           % ELSEIF the coefficient is Beta
           elseif (strcmp(Continuation.Coefficient_name{m},'Beta'))
               Physics2D = Physics2D_Var2d(Method,Physics2D.Delta,Continuation.Coefficient{m}{n},Physics2D.Omega); % Replace the coefficient Beta
           % ELSEIF the coefficient is Omega
           elseif (strcmp(Continuation.Coefficient_name{m},'Omega'))
               Physics2D = Physics2D_Var2d(Method,Physics2D.Delta,Physics2D.Beta,Continuation.Coefficient{m}{n}); % Replace the coefficient Omega
           % ELSEIF the coefficient is the G matrix of the potential
           elseif (strcmp(Continuation.Coefficient_name{m},'GPotential'))
               Physics2D = Potential_Var2d(Method, Physics2D, TMP_Potential_function, Continuation.Coefficient{m}{n}); % Replace the G matrix of the potential
           % ELSEIF the coefficient is the G matrix of the nonlinearity
           elseif (strcmp(Continuation.Coefficient_name{m},'GNonlinearity'))
               Physics2D = Nonlinearity_Var2d(Method, Physics2D, TMP_Nonlinearity_function, Continuation.Coefficient{m}{n}, TMP_Nonlinearity_energy_function); % Replace the G matrix of the nonlinearity
           % ELSEIF the coefficient is the G matrix of the non-local nonlinearity
           elseif (strcmp(Continuation.Coefficient_name{m},'GFFTNonlinearity'))
               Physics2D = FFTNonlinearity_Var2d(Method, Physics2D, TMP_FFTNonlinearity_function, Continuation.Coefficient{m}{n}, TMP_FFTNonlinearity_energy_function); % Replace the G matrix of the non-local nonlinearity
           % ELSEIF the coefficient is the G matrix of the gradientx
           elseif (strcmp(Continuation.Coefficient_name{m},'GGradientx'))
               Physics2D = Gradientx_Var2d(Method, Physics2D, TMP_Gradientx_function, Continuation.Coefficient{m}{n}); % Replace the G matrix of the gradientx
           % ELSEIF the coefficient is the G matrix of the gradienty
           elseif (strcmp(Continuation.Coefficient_name{m},'GGradienty'))
               Physics2D = Gradienty_Var2d(Method, Physics2D, TMP_Gradienty_function, Continuation.Coefficient{m}{n}); % Replace the G matrix of the gradienty
           end
       end
    end
    %% Ground state computation
    if (strcmp(Method.Computation,'Ground'))
        if (strcmp(Method.Type,'BESP'))
        [Phi,Outputs] = BESP_CNGF2d(Phi_0, Method, Geometry2D, Physics2D, Outputs, Print, Figure);
        elseif (strcmp(Method.Type,'AccBESP'))
        [Phi,Outputs] = AccBESP_CNGF2d(Phi_0, Method, Geometry2D, Physics2D, Outputs, Print, Figure);
        elseif (strcmp(Method.Type,'CNSP'))
        [Phi,Outputs] = CNSP_CNGF2d(Phi_0, Method, Geometry2D, Physics2D, Outputs, Print, Figure);
        elseif (strcmp(Method.Type,'BEFD'))
        [Phi, Outputs] = BEFD_CNGF2d(Phi_0, Method, Geometry2D, Physics2D, Outputs, Print, Figure);
        elseif (strcmp(Method.Type,'CNFD'))
        [Phi, Outputs] = CNFD_CNGF2d(Phi_0, Method, Geometry2D, Physics2D, Outputs, Print, Figure);
        end
    %% Dynamic solution computation
    elseif (strcmp(Method.Computation,'Dynamic'))
        if (strcmp(Method.Type,'BESP'))
        [Phi,Outputs] = BESP_DS2d(Phi_0, Method, Geometry2D, Physics2D, Outputs, Print, Figure);
        elseif (strcmp(Method.Type,'CNSP'))
        [Phi,Outputs] = CNSP_DS2d(Phi_0, Method, Geometry2D, Physics2D, Outputs, Print, Figure);
        elseif (strcmp(Method.Type,'BEFD'))
        [Phi, Outputs] = BEFD_DS2d(Phi_0, Method, Geometry2D, Physics2D, Outputs, Print, Figure);
        elseif (strcmp(Method.Type,'CNFD'))
        [Phi, Outputs] = CNFD_DS2d(Phi_0, Method, Geometry2D, Physics2D, Outputs, Print, Figure);
        elseif (strcmp(Method.Type,'Relaxation'))
        [Phi, Outputs] = RSP_DS2d(Phi_0, Method, Geometry2D, Physics2D, Outputs, Print, Figure);
        elseif (strcmp(Method.Type,'Splitting'))
        [Phi, Outputs] = SPL_DS2d(Phi_0, Method, Geometry2D, Physics2D, Outputs, Print, Figure);
        end
    end

end