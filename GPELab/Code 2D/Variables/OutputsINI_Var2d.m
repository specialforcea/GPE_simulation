%% Initialization of the Outputs variables
%% INPUT:
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%% INPUTS(OPTIONAL):
%%          Evo_outputs: Variable containing the number of iterations between each computation of the outputs
%%          save: Variable containing the information wherever to save or not the solution during the computation
%%          outputs: User defined functions for each component which depend on (phi,X,Y) to be stored and printed in the informations (cell array of functions)
%%          outputs_names: Names of the user defined functions (cell array of char)
%% OUTPUT:
%%          Outputs: Different outputs computed during the computation of the ground states (structure)

function [Outputs] = OutputsINI_Var2d(Method, Evo_outputs, save, userdef_outputs, userdef_outputs_names, globaluserdef_outputs, globaluserdef_outputs_names)
%% Initialization of the Outputs variables
    Outputs.phi_abs_0 = cell(1,Method.Ncomponents); % Initialization of the variable containing the square of the wave function at the origin
    Outputs.x_rms = cell(1,Method.Ncomponents); % Initialization of the variable containing the root mean square in the x direction
    Outputs.y_rms = cell(1,Method.Ncomponents); % Initialization of the variable containing the root mean square in the y direction
    Outputs.Energy = cell(1,Method.Ncomponents); % Initialization of the variable containing the energy of the wave function 
    Outputs.Chemical_potential = cell(1,Method.Ncomponents); % Initialization of the variable containing the chemical potential of the wave function 
    Outputs.Angular_momentum = cell(1,Method.Ncomponents); % Initialization of the variable containing the angular momentum of the wave function
    Outputs.Iterations = 0; % Initialization of the variable containing the index of calculated outputs
    Outputs.Save_solution = 0; % Initialization of the variable containing the information wherever to save or no the solution during the computation
    % If one has set the evo variable
    if (nargin >=2) && (isempty(Evo_outputs) == 0)
        % IF one has chosen to save the functions during the simulation
        if (isposintscalar(Evo_outputs))
            Outputs.Evo_outputs = Evo_outputs; % Initialization of the variable containing the number of iterations between each computation of the outputs
        else
            Outputs.Evo_outputs = 5; % Initialization of the variable containing the number of iterations between each computation of the outputs
        end
    else
         Outputs.Evo_outputs = 5; % Initialization of the variable containing the number of iterations between each computation of the outputs
    end

    % IF one wants to save the wave functions during the simulation
    if (nargin >=3) && (isempty(save) == 0)
        % IF one has chosen to save the functions during the simulation
        if (save == 1)
            Outputs.Save_solution = 1; % Initialization of the variable containing the information wherever to save or not the solution during the simulation
        else
            Outputs.Save_solution = 0; % Initialization of the variable containing the information wherever to save or not the solution during the simulation
        end
    else
        Outputs.Save_solution = 0; % Initialization of the variable containing the information wherever to save or not the solution during the simulation
    end
% IF there are user defined outputs local functions
if (nargin >= 4) && (isempty(userdef_outputs) == 0)
%% Initialization of the user defined outputs
    Outputs.User_defined_number_local = length(userdef_outputs); % Initialization of the number of user defined functions
    Outputs.User_defined_function_local = userdef_outputs; % Initialization of the variable containing the user defined functions
    Outputs.User_compute_local = 1; % Setting to compute user defined functions
    % IF there are names for the user defined outputs functions
    if (nargin >= 5) && (isempty(userdef_outputs_names) == 0)
        Outputs.User_defined_names_local = userdef_outputs_names; % Storing the names of the user defined functions
    % ELSE if there are no names for the user defined outputs functions
    else
        % FOR each user defined function
        for n = 1:Outputs.User_defined_number_local
            Outputs.User_defined_names_local{n} = strcat('User defined local function',32,num2str(n)); % Storing the default name
        end
    end
% ELSE if there are no user defined outputs functions
else
    Outputs.User_compute_local = 0; % Setting not to compute user defined local functions
end
    % IF there are user defined outputs global functions
if (nargin >= 6) && (isempty(globaluserdef_outputs) == 0)
%% Initialization of the user defined outputs
    Outputs.User_defined_number_global = length(globaluserdef_outputs); % Initialization of the number of user defined functions
    Outputs.User_defined_function_global = globaluserdef_outputs; % Initialization of the variable containing the user defined functions
    Outputs.User_compute_global = 1; % Setting to compute user defined functions
    % IF there are names for the user defined outputs functions
    if (nargin == 7) && (isempty(globaluserdef_outputs_names) == 0)
        Outputs.User_defined_names_global = globaluserdef_outputs_names; % Storing the names of the user defined functions
    % ELSE if there are no names for the user defined outputs functions
    else
        % FOR each user defined function
        for n = 1:Outputs.User_defined_number_global
            Outputs.User_defined_names_global{n} = strcat('User defined global function',32,num2str(n)); % Storing the default name
        end
    end
% ELSE if there are no user defined outputs functions
else
    Outputs.User_compute_global = 0; % Setting not to compute user defined global functions
end