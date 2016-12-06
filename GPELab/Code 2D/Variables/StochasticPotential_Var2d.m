%% Addition of the stochastic potential function matrix in the physics of the
%% problem
%% INPUTS:
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          Physics2D: Structure containing variables concerning the physics of the problem in 2D (structure) (see Physics2D_Var2d.m)
%% INPUTS(OPTIONAL):
%%          Potential: Function or cell array of functions that depends on (W,X,Y) which will be multiplied by the wave function in the physical problem (function or cell array of functions)
%%          (In the case of a function, the function will be applied only on the diagonal terms)
%%          G: Matrix that will be multiplied element by element to the potential function matrix (matrix)
%%          StochasticProcess: Function that depends on (t) which gives the increment between time t-dt and t and which will be used as the stochastic process for the simulation (function)
%% OUTPUT:
%%          Physics2D: Structure containing variables concerning the physics of the problem in 2D (structure) (see Physics2D_Var2d.m)

function [Physics2D] = StochasticPotential_Var2d(Method, Physics2D, StochasticPotential, G, StochasticProcess)
%% Initializing the default potential
Default_Potential = @(W,X,Y)quadratic_potential2d(1,1,X,Y);

%% Adding the stochastic process to the Physics2D variable
[Process,Increments] = Brownian_Process2d(Method,1/2);
if (nargin == 5) && (isempty(StochasticProcess) == 0) && (iscell(StochasticProcess) == 0)
    Physics2D.StochasticProcess_function = @(t,X,Y) StochasticProcess(t,X,Y);
elseif (nargin == 5) && (isempty(StochasticProcess) == 0) && (iscell(StochasticProcess) == 1)
    for m = 1:length(StochasticProcess)
        Physics2D.StochasticProcess_function{m} = @(t,X,Y) StochasticProcess{m}(t,X,Y);
    end
else
    Physics2D.StochasticProcess_function = @(t,X,Y) Process(t);
end

%% Adding the potential function matrix
% IF there are 5 inputs and the potential is not defined
if (nargin == 5) && (iscell(StochasticPotential) == 0) && (isempty(StochasticPotential) == 1) && (isempty(G) == 0)
    % FOR each component
    for n = 1:Method.Ncomponents
        Potential_function_index = []; % Initializing the temporary potential index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero
                if (G(n,m)~= 0)
                    % IF it is a diagonal term
                    if(n == m)
                    Physics2D.StochasticPotential_function{n,m} = @(W,X,Y)G(n,m)*Default_Potential(W,X,Y); % Storing the potential as the correspondant value of G at the index added by the default potential
                    Potential_function_index = [Potential_function_index,m]; % Add the 'm' index in the temporary potential index
                    elseif (n ~= m)
                    Physics2D.StochasticPotential_function{n,m} = @(W,X,Y) 0; % Storing zero in the potential function matrix
                    end
                % ELSEIF the value of the matrix at the index is zero
                elseif (G(n,m) == 0)
                    Physics2D.StochasticPotential_function{n,m} = @(W,X,Y) 0; % Storing zero in the potential function matrix
                end
        end
        Physics2D.StochasticPotential_function_Index{n} = Potential_function_index; % Store the potential index for the 'm' index
    end
% ELSEIF there are 5 inputs and the potential is defined but not a
% cell array
elseif (nargin == 5) && (iscell(StochasticPotential) == 0) && (isempty(StochasticPotential) == 0) && (isempty(G) == 0)
    % FOR each component
    for n = 1:Method.Ncomponents
        Potential_function_index = []; % Initializing the temporary potential index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero
                if (G(n,m)~= 0)
                    Physics2D.StochasticPotential_function{n,m} = @(W,X,Y)G(n,m)*StochasticPotential(W,X,Y); % Storing the potential as the correspondant value of G at the index multiplied by the defined potential
                    Potential_function_index = [Potential_function_index,m]; % Add the 'm' index in the temporary potential index
                % ELSEIF the value of the matrix at the index is zero
                elseif (G(n,m) == 0)
                    Physics2D.StochasticPotential_function{n,m} = @(W,X,Y)0; % Storing zero in the potential function matrix
                end
        end
        Physics2D.StochasticPotential_function_Index{n} = Potential_function_index; % Store the potential index for the 'm' index
    end
% ELSEIF there are 5 inputs and the potential is a cell array
elseif (nargin == 5) && (iscell(StochasticPotential) == 1) && (isempty(G) == 0)
    % FOR each component
    for n = 1:Method.Ncomponents
        Potential_function_index = []; % Initializing the temporary potential index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero and
                % the potential function is defined
                if (G(n,m)~= 0) && (isempty(StochasticPotential{n,m}) == 0)
                    Physics2D.StochasticPotential_function{n,m} = @(W,X,Y)G(n,m)*StochasticPotential{n,m}(W,X,Y); % Storing the potential as the correspondant value of G at the index multiplied by the defined potential
                    Potential_function_index = [Potential_function_index,m]; % Add the 'm' index in the temporary potential index
                % ELSEIF the value of the matrix at the index is zero or
                % the potential function is not defined
                else
                    Physics2D.StochasticPotential_function{n,m} = @(W,X,Y)0; % Storing zero in the potential function matrix
                end
        end
        Physics2D.StochasticPotential_function_Index{n} = Potential_function_index; % Store the potential index for the 'm' index
    end
    
% ELSEIF there are 5 inputs and the potential is not defined
elseif (nargin == 5) && (iscell(StochasticPotential) == 0) && (isempty(StochasticPotential) == 1) && (isempty(G) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Potential_function_index = []; % Initializing the temporary potential index
        % FOR each component
        for m = 1:Method.Ncomponents
                    % IF it is a diagonal term
                    if(n == m)
                    Physics2D.StochasticPotential_function{n,m} = @(W,X,Y)Default_Potential(W,X,Y); % Storing the potential as the correspondant value of G at the index added by the default potential
                    Potential_function_index = [Potential_function_index,m]; % Add the 'm' index in the temporary potential index
                    % ELSE if it is an extradiagonal term
                    else
                    Physics2D.StochasticPotential_function{n,m} = @(W,X,Y) 0; % Storing zero in the potential function matrix
                    end
        end
        Physics2D.StochasticPotential_function_Index{n} = Potential_function_index; % Store the potential index for the 'm' index
    end
% ELSEIF there are 5 inputs and the potential is defined but not a
% cell array
elseif (nargin == 5) && (iscell(StochasticPotential) == 0) && (isempty(StochasticPotential) == 0) && (isempty(G) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Potential_function_index = []; % Initializing the temporary potential index
        % FOR each component
        for m = 1:Method.Ncomponents
                    % IF it is a diagonal term
                    if (n == m)
                    Physics2D.StochasticPotential_function{n,m} = @(W,X,Y)StochasticPotential(W,X,Y); % Storing the potential as the correspondant value of G at the index multiplied by the defined potential
                    Potential_function_index = [Potential_function_index,m]; % Add the 'm' index in the temporary potential index
                    % ELSE if it is an extradiagonal term
                    else
                    Physics2D.StochasticPotential_function{n,m} = @(W,X,Y)0; % Storing zero in the potential function matrix
                    end
        end
        Physics2D.StochasticPotential_function_Index{n} = Potential_function_index; % Store the potential index for the 'm' index
    end
% ELSEIF there are 5 inputs and the potential is a cell array
elseif (nargin == 5) && (iscell(StochasticPotential) == 1) && (isempty(G) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Potential_function_index = []; % Initializing the temporary potential index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero and
                % the potential function is defined
                if (isempty(StochasticPotential{n,m}) == 0)
                    Physics2D.StochasticPotential_function{n,m} = @(W,X,Y)StochasticPotential{n,m}(W,X,Y); % Storing the potential as the correspondant value of G at the index multiplied by the defined potential
                    Potential_function_index = [Potential_function_index,m]; % Add the 'm' index in the temporary potential index
                % ELSEIF the value of the matrix at the index is zero or
                % the potential function is not defined
                else
                    Physics2D.StochasticPotential_function{n,m} = @(W,X,Y)0; % Storing zero in the potential function matrix
                end
        end
        Physics2D.StochasticPotential_function_Index{n} = Potential_function_index; % Store the potential index for the 'm' index
    end
% ELSEIF there are 4 inputs and the potential is not defined
elseif (nargin == 4) && (iscell(StochasticPotential) == 0) && (isempty(StochasticPotential) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Potential_function_index = []; % Initializing the temporary potential index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero
                if (G(n,m)~= 0)
                    % IF it is a diagonal term
                    if(n == m)
                    Physics2D.StochasticPotential_function{n,m} = @(W,X,Y)G(n,m)*Default_Potential(W,X,Y); % Storing the potential as the correspondant value of G at the index added by the default potential
                    Potential_function_index = [Potential_function_index,m]; % Add the 'm' index in the temporary potential index
                    elseif (n ~= m)
                    Physics2D.StochasticPotential_function{n,m} = @(W,X,Y) 0; % Storing zero in the potential function matrix
                    end
                % ELSEIF the value of the matrix at the index is zero
                elseif (G(n,m) == 0)
                    Physics2D.StochasticPotential_function{n,m} = @(W,X,Y) 0; % Storing zero in the potential function matrix
                end
        end
        Physics2D.StochasticPotential_function_Index{n} = Potential_function_index; % Store the potential index for the 'm' index
    end
% ELSEIF there are 4 inputs and the potential is defined but not a
% cell array
elseif (nargin == 4) && (iscell(StochasticPotential) == 0) && (isempty(StochasticPotential) == 0)
    % FOR each component
    for n = 1:Method.Ncomponents
        Potential_function_index = []; % Initializing the temporary potential index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero
                if (G(n,m)~= 0)
                    Physics2D.StochasticPotential_function{n,m} = @(W,X,Y)G(n,m)*StochasticPotential(W,X,Y); % Storing the potential as the correspondant value of G at the index multiplied by the defined potential
                    Potential_function_index = [Potential_function_index,m]; % Add the 'm' index in the temporary potential index
                % ELSEIF the value of the matrix at the index is zero
                elseif (G(n,m) == 0)
                    Physics2D.StochasticPotential_function{n,m} = @(W,X,Y)0; % Storing zero in the potential function matrix
                end
        end
        Physics2D.StochasticPotential_function_Index{n} = Potential_function_index; % Store the potential index for the 'm' index
    end
% ELSEIF there are 4 inputs and the potential is a cell array
elseif (nargin == 4) && (iscell(StochasticPotential) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Potential_function_index = []; % Initializing the temporary potential index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero and
                % the potential function is defined
                if (G(n,m)~= 0) && (isempty(StochasticPotential{n,m}) == 0)
                    Physics2D.StochasticPotential_function{n,m} = @(W,X,Y)G(n,m)*StochasticPotential{n,m}(W,X,Y); % Storing the potential as the correspondant value of G at the index multiplied by the defined potential
                    Potential_function_index = [Potential_function_index,m]; % Add the 'm' index in the temporary potential index
                % ELSEIF the value of the matrix at the index is zero or
                % the potential function is not defined
                else
                    Physics2D.StochasticPotential_function{n,m} = @(W,X,Y)0; % Storing zero in the potential function matrix
                end
        end
        Physics2D.StochasticPotential_function_Index{n} = Potential_function_index; % Store the potential index for the 'm' index
    end
% ELSEIF there are 3 inputs and the potential is not defined
elseif (nargin == 3) && (iscell(StochasticPotential) == 0) && (isempty(StochasticPotential) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Potential_function_index = []; % Initializing the temporary potential index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF it is a diagonal term
                if (n == m)
                    Physics2D.StochasticPotential_function{n,m} = @(W,X,Y)Default_Potential(W,X,Y); % Storing the potential as the default potential
                    Potential_function_index = [Potential_function_index,m]; % Add the 'm' index in the temporary potential index
                % ELSEIF it is an extradiagonal term
                elseif (n ~= m)
                    Physics2D.StochasticPotential_function{n,m} = @(W,X,Y)0; % Storing zero in the potential function matrix
                end
        end
        Physics2D.StochasticPotential_function_Index{n} = Potential_function_index; % Store the potential index for the 'm' index
    end
% ELSEIF there are 3 inputs and the potential is defined but not a
% cell array
elseif (nargin == 3) && (iscell(StochasticPotential) == 0) && (isempty(StochasticPotential) == 0)
    % FOR each component
    for n = 1:Method.Ncomponents
        Potential_function_index = []; % Initializing the temporary potential index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF it is a diagonal term
                if (n == m)
                    Physics2D.StochasticPotential_function{n,m} = @(W,X,Y)StochasticPotential(W,X,Y); % Storing the potential as the defined potential
                    Potential_function_index = [Potential_function_index,m]; % Add the 'm' index in the temporary potential index
                % ELSEIF it is an extradiagonal term
                elseif (n ~= m)
                    Physics2D.StochasticPotential_function{n,m} = @(W,X,Y)0; % Storing zero in the potential function matrix
                end
        end
        Physics2D.StochasticPotential_function_Index{n} = Potential_function_index; % Store the potential index for the 'm' index
    end
% ELSEIF there are 3 inputs and the potential is a cell array
elseif (nargin == 3) && (iscell(StochasticPotential) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Potential_function_index = []; % Initializing the temporary potential index
        % FOR each component
        for m = 1:Method.Ncomponents
            % IF the gradient function is defined
            if (isempty(StochasticPotential{n,m}) == 0)
                Physics2D.StochasticPotential_function{n,m} = @(W,X,Y)StochasticPotential{n,m}(W,X,Y); % Storing the potential as the defined potential
                Potential_function_index = [Potential_function_index,m]; % Add the 'm' index in the temporary potential index
            %ELSE if the gradient function is not defined
            else
                Physics2D.StochasticPotential_function{n,m} = @(W,X,Y) 0; % Storing zero in the potential function matrix
            end
        end
        Physics2D.StochasticPotential_function_Index{n} = Potential_function_index; % Store the potential index for the 'm' index
    end
% ELSEIF there are 2 inputs
elseif (nargin == 2)
    % FOR each component
    for n = 1:Method.Ncomponents
        Potential_function_index = []; % Initializing the temporary potential index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF it is a diagonal term
                if (n == m)
                    Physics2D.StochasticPotential_function{n,m} = @(W,X,Y)Default_Potential(W,X,Y); % Storing the potential as the default potential
                    Potential_function_index = [Potential_function_index,m]; % Add the 'm' index in the temporary potential index
                % ELSEIF it is an extradiagonal term
                elseif (n ~= m)
                    Physics2D.StochasticPotential_function{n,m} = @(W,X,Y)0; % Storing zero in the potential function matrix
                end
        end
        Physics2D.StochasticPotential_function_Index{n} = Potential_function_index; % Store the potential index for the 'm' index
    end
end