%% Launching simulation and applying the continuation method
%% INPUTS:
%%          Phi_0: Initial wave functions (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var1d.m)
%%          Geometry1D: Structure containing variables concerning the geometry of the problem in 1D (structure) (see Geometry1D_Var1d.m)
%%          Physics1D: Structure containing variables concerning the physics of the problem in 1D (structure) (see Physics1D_Var1d.m)
%%          Outputs: Different outputs computed during the computation of the ground states (structure) (see OutputsINI_Var1D.m and FFTOutputs_Var1d.m)
%% INPUTS(OPTIONAL):
%%          Continuation: Continuation structure used to apply a continuation method (structure) (see Continuation_Var1d.m)
%%          Print: Structure containing variables concerning the printing and drawing of informations during the computations (structure) (see Print_Var1d.m)
%%          Figure: Structure containing variables concerning the figures (structure)
%% OUTPUTS:
%%          Phi: Wave function at the end of the simulation (cell array)
%%          Outputs: Different outputs computed during the simulation (structure) (see OutputsINI_Var1d.m and FFTOutputs_Var1d.m)
%% FUNCTIONS USED:
%%          Physics1D_Var1d: To replace the Delta, Beta or Omega coefficients for the continuation method (line 85 and 88)
%%          Potential_Var1d: To replace the G matrix of the potential for the continuation method (line 91)
%%          Nonlinearity_Var1d: To replace the G matrix of the nonlinearity for the continuation method (line 94)
%%          Gradientx_Var1d: To replace the G matrix of the gradientx for the continuation method (line 97)
%%          BESP_CNGF1d: To compute a ground state using the BESP scheme (line 104)
%%          CNSP_CNGF1d: To compute a ground state using the CNSP scheme (line 106)
%%          BEFD_CNGF1d: To compute a ground state using the BEFD scheme (line 108)
%%          CNFD_CNGF1d: To compute a ground state using the CNFD scheme (line 110)
%%          BESP_DS1d: To compute a dynamic solution using the BESP scheme (line 115)
%%          CNSP_DS1d: To compute a dynamic solution using the CNSP scheme (line 117)
%%          BEFD_DS1d: To compute a dynamic solution using the BEFD scheme (line 119)
%%          CNFD_DS1d: To compute a dynamic solution using the CNFD scheme (line 121)

function [Phi, Outputs] = GPELab1d(varargin)
%% Analysis of the inputs
Analyse_Var = inputParser; % Creating the parser
Analyse_Var.addRequired('Phi_0',@(x)iscell(x)); % Required input 'Phi_0' which must be a matrix
Analyse_Var.addRequired('Method',@(x)isstruct(x)); % Required input 'Method' which must be a structure
Analyse_Var.addRequired('Geometry1D',@(x)isstruct(x)); % Required input 'Geometry1D' which must be a structure
Analyse_Var.addRequired('Physics1D',@(x)isstruct(x)); % Required input 'Physics1D' which must be a structure
Analyse_Var.addRequired('Outputs',@(x)isstruct(x)); % Required input 'Outputs' which must be a structure
Analyse_Var.addOptional('Continuation',[],@(x)isstruct(x)+isempty(x)); % Optional input 'Continuation' which must be a structure
Analyse_Var.addOptional('Print', Print_Var1d(1,10,0),@(x)isstruct(x)); % Optional input 'Print' which must be a structure
Analyse_Var.addOptional('Figure', Figure_Var1d,@(x)isstruct(x)); % Optional input 'Figure' which must be a structure

%% Parsing inputs and creating the Method structure
% Parsing inputs
Analyse_Var.parse(varargin{:}); % Analysing the inputs
% Contructing the Method structure
Phi_0 = Analyse_Var.Results.Phi_0; %Storing the 'Phi_0' input
Method = Analyse_Var.Results.Method; %Storing the 'Method' input
Geometry1D = Analyse_Var.Results.Geometry1D; %Storing the 'Geometry1D' input
Physics1D = Analyse_Var.Results.Physics1D; %Storing the 'Physics1D' input
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
TMP_Potential_function = Physics1D.Potential_function; % Storing temporary potential functions
TMP_Nonlinearity_function = Physics1D.Nonlinearity_function; % Storing temporary nonlinearity functions
TMP_Nonlinearity_energy_function = Physics1D.Nonlinearity_energy_function; % Storing temporary nonlinearity functions
TMP_FFTNonlinearity_function = Physics1D.FFTNonlinearity_function; % Storing temporary non-local nonlinearity functions
TMP_FFTNonlinearity_energy_function = Physics1D.FFTNonlinearity_energy_function; % Storing temporary non-local nonlinearity functions
TMP_Gradientx_function = Physics1D.Gradientx_function; % Storing temporary gradientx functions
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
               Physics1D = Physics1D_Var1d(Method, Continuation.Coefficient{m}{n},Physics1D.Beta); % Replace the coefficient Delta
           % ELSEIF the coefficient is Beta
           elseif (strcmp(Continuation.Coefficient_name{m},'Beta'))
               Physics1D = Physics1D_Var1d(Method,Physics1D.Delta,Continuation.Coefficient{m}{n}); % Replace the coefficient Beta
           % ELSEIF the coefficient is the G matrix of the potential
           elseif (strcmp(Continuation.Coefficient_name{m},'GPotential'))
               Physics1D = Potential_Var1d(Method, Physics1D, TMP_Potential_function, Continuation.Coefficient{m}{n}); % Replace the G matrix of the potential
           % ELSEIF the coefficient is the G matrix of the nonlinearity
           elseif (strcmp(Continuation.Coefficient_name{m},'GNonlinearity'))
               Physics1D = Nonlinearity_Var1d(Method, Physics1D, TMP_Nonlinearity_function, Continuation.Coefficient{m}{n}, TMP_Nonlinearity_energy_function); % Replace the G matrix of the nonlinearity
           % ELSEIF the coefficient is the G matrix of the non-local nonlinearity
           elseif (strcmp(Continuation.Coefficient_name{m},'GFFTNonlinearity'))
               Physics1D = FFTNonlinearity_Var1d(Method, Physics1D, TMP_FFTNonlinearity_function, Continuation.Coefficient{m}{n}, TMP_FFTNonlinearity_energy_function); % Replace the G matrix of the non-local nonlinearity
           % ELSEIF the coefficient is the G matrix of the gradientx
           elseif (strcmp(Continuation.Coefficient_name{m},'GGradientx'))
               Physics1D = Gradientx_Var1d(Method, Physics1D, TMP_Gradientx_function, Continuation.Coefficient{m}{n}); % Replace the G matrix of the gradientx
           end
       end
    end
    %% Ground state computation
    if (strcmp(Method.Computation,'Ground'))
        if (strcmp(Method.Type,'BESP'))
        [Phi,Outputs] = BESP_CNGF1d(Phi_0, Method, Geometry1D, Physics1D, Outputs, Print, Figure);
        elseif (strcmp(Method.Type,'CNSP'))
        [Phi,Outputs] = CNSP_CNGF1d(Phi_0, Method, Geometry1D, Physics1D, Outputs, Print, Figure);
        elseif (strcmp(Method.Type,'BEFD'))
        [Phi, Outputs] = BEFD_CNGF1d(Phi_0, Method, Geometry1D, Physics1D, Outputs, Print, Figure);
        elseif (strcmp(Method.Type,'CNFD'))
        [Phi, Outputs] = CNFD_CNGF1d(Phi_0, Method, Geometry1D, Physics1D, Outputs, Print, Figure);
        end
    %% Dynamic solution computation
    elseif (strcmp(Method.Computation,'Dynamic'))
        if (strcmp(Method.Type,'BESP'))
        [Phi,Outputs] = BESP_DS1d(Phi_0, Method, Geometry1D, Physics1D, Outputs, Print, Figure);
        elseif (strcmp(Method.Type,'CNSP'))
        [Phi,Outputs] = CNSP_DS1d(Phi_0, Method, Geometry1D, Physics1D, Outputs, Print, Figure);
        elseif (strcmp(Method.Type,'BEFD'))
        [Phi, Outputs] = BEFD_DS1d(Phi_0, Method, Geometry1D, Physics1D, Outputs, Print, Figure);
        elseif (strcmp(Method.Type,'CNFD'))
        [Phi, Outputs] = CNFD_DS1d(Phi_0, Method, Geometry1D, Physics1D, Outputs, Print, Figure);
        elseif (strcmp(Method.Type,'Relaxation'))
        [Phi, Outputs] = RSP_DS1d(Phi_0, Method, Geometry1D, Physics1D, Outputs, Print, Figure);
        elseif (strcmp(Method.Type,'Splitting'))
        [Phi, Outputs] = SPL_DS1d(Phi_0, Method, Geometry1D, Physics1D, Outputs, Print, Figure);
        end
    end

end