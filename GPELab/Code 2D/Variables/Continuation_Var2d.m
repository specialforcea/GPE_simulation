%% Creation of the continuation structure
%% INPUTS:
%%          Coefficient: Evolution of the coefficient during each computation of the continuation (cell vector of matrices or cell vector of doubles)
%%          Coefficient_name: Name associated to the coefficient (string)
%%          (Must either be: 'Delta','Beta','Omega','GPotential','GNonlinearity','GGradientx','GGradienty')
%% INPUTS(OPTIONAL):
%%          Continuation: Continuation structure (structure, Default: [])
%% OUTPUT:
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)

function [Continuation] = Continuation_Var2d(varargin)
%% Analysis of the inputs
valid_Coefficient_name = {'Delta','Beta','Omega','GPotential','GNonlinearity','GFFTNonlinearity','GGradientx','GGradienty'}; % List of valid inputs for the coefficient name
Analyse_Var = inputParser; % Creating the parser
Analyse_Var.addRequired('Coefficient_name',@(x)ischar(validatestring(x,valid_Coefficient_name))); % Required input 'Coefficient_name' which must be a valid string
Analyse_Var.addRequired('Coefficient',@(x)iscell(x)); % Required input 'Coefficient' which must be numeric
Analyse_Var.addOptional('Continuation',[],@(x)isstruct(x)); % Optional input 'Continuation' which must be a structure

%% Parsing inputs
% Parsing inputs
Analyse_Var.parse(varargin{:}); % Analysing the inputs
% Storing the inputs
Coefficient_name = Analyse_Var.Results.Coefficient_name; %Storing the 'Coefficient_name' input
Coefficient = Analyse_Var.Results.Coefficient; %Storing the 'Coefficient' input
Continuation_tmp = Analyse_Var.Results.Continuation; %Storing the 'Continuation' input

%% Building or updating the Continuation structure
% IF the continuation structure is not given
if (isempty(Continuation_tmp))
    Continuation.Ncontinuation = 1; % Initializing the number of coefficients
    Continuation.Coefficient{Continuation.Ncontinuation} = Coefficient; % Storing the coefficient
    Continuation.Coefficient_name{Continuation.Ncontinuation} = Coefficient_name; % Storing the coefficient's name
% ELSE if the continuation structure is not given
else
    Continuation = Continuation_tmp; % Storing the continuation structure
    Rewrite_num = 0; % Initializing a coordinate in case of rewrite
    % FOR each coefficient name
    for n = 1:length(Continuation.Coefficient_name)
        % IF the name of the coefficient is the same as the one we want to
        % add
        if (strcmp(Continuation.Coefficient_name{n},Coefficient_name))
            Rewrite_num = n; % Storing the coordinate of the coefficient
        end
    end
    % IF there is no need to rewrite the coefficient
    if (Rewrite_num == 0)
        Continuation.Ncontinuation = Continuation.Ncontinuation + 1; % Updating the number of coefficients
        Continuation.Coefficient{Continuation.Ncontinuation} = Coefficient; % Storing the coefficient
        Continuation.Coefficient_name{Continuation.Ncontinuation} = Coefficient_name; % Storign the coefficient's name
    % ELSE if the coefficient is the same, we rewrite the coefficient
    else
        Continuation.Coefficient{Rewrite_num} = Coefficient; % Rewriting the coefficient
        Continuation.Coefficient_name{Rewrite_num} = Coefficient_name; % Rewritinh the coefficient's name
    end
end