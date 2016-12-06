%% Addition of the gradient in the y direction function matrix in the physics of the
%% problem
%% INPUTS:
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          Physics2D: Structure containing variables concerning the physics of the problem in 2D (structure) (see Physics2D_Var2d.m)
%% INPUTS(OPTIONAL):
%%          Gradienty: Function or cell array of functions that depends on (X,Y) which will be multiplied by the gradient in the y direction of the wave function in the physical problem (function or cell array of functions)
%%          (In the case of a function, the function will be applied only on the diagonal terms)
%%          G: Matrix that will be multiplied element by element to the gradient in the y direction function matrix (matrix)
%% OUTPUT:
%%          Physics2D: Structure containing variables concerning the physics of the problem in 2D (structure) (see Physics2D_Var2d.m)

function [Physics2D] = Gradienty_Var2d(Method, Physics2D, Gradienty, G)
%% Initializing the default gradient function
Default_Gradienty = @(X,Y) -1i*Physics2D.Omega*X; % Setting the default gradient function

%% Adding the gradient in the y direction function matrix
% IF there are 4 inputs and the gradient function is not defined
if (nargin == 4) && (iscell(Gradienty) == 0) && (isempty(Gradienty) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Gradienty_function_index = []; % Initializing the temporary gradient index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero
                if (G(n,m)~= 0)
                    % IF it is a diagonal term
                    if(n == m)
                    Physics2D.Gradienty_function{n,m} = @(X,Y)G(n,m)*Default_Gradienty(X,Y); % Storing the gradient function as the correspondant value of G at the index added by the default gradient function
                    Gradienty_function_index = [Gradienty_function_index,m]; % Add the 'm' index in the temporary gradient function index
                    Physics2D.Gradienty_compute_Index{m} = 1; % Setting the computation gradient index
                    % ELSE if it is an extradiagonal term
                    else
                    Physics2D.Gradienty_function{n,m} = @(X,Y) 0; % Storing zero in the gradient function matrix
                    end
                % ELSEIF the value of the matrix at the index is zero
                elseif (G(n,m) == 0)
                    Physics2D.Gradienty_function{n,m} = @(X,Y) 0; % Storing zero in the gradient function matrix
                end
        end
        Physics2D.Gradienty_function_Index{n} = Gradienty_function_index; % Store the gradient function index for the 'm' index
    end
% ELSEIF there are 4 inputs and the gradient function is defined but not a
% cell array
elseif (nargin == 4) && (iscell(Gradienty) == 0) && (isempty(Gradienty) == 0)
    % FOR each component
    for n = 1:Method.Ncomponents
        Gradienty_function_index = []; % Initializing the temporary gradient function index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero
                if (G(n,m)~= 0)
                    Physics2D.Gradienty_function{n,m} = @(X,Y)G(n,m)*Gradienty(X,Y); % Storing the gradient function as the correspondant value of G at the index multiplied by the defined gradient function
                    Gradienty_function_index = [Gradienty_function_index,m]; % Add the 'm' index in the temporary gradient function index
                    Physics2D.Gradienty_compute_Index{m} = 1; % Setting the computation gradient index
                % ELSEIF the value of the matrix at the index is zero
                elseif (G(n,m) == 0)
                    Physics2D.Gradienty_function{n,m} = @(X,Y)0; % Storing zero in the gradient function matrix
                end
        end
        Physics2D.Gradienty_function_Index{n} = Gradienty_function_index; % Store the gradient function index for the 'm' index
    end
% ELSEIF there are 4 inputs and the gradient function is a cell array
elseif (nargin == 4) && (iscell(Gradienty) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Gradienty_function_index = []; % Initializing the temporary gradient function index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero and
                % the gradient function is defined
                if (G(n,m)~= 0) && (isempty(Gradienty{n,m}) == 0)
                    Physics2D.Gradienty_function{n,m} = @(X,Y)G(n,m)*Gradienty{n,m}(X,Y); % Storing the gradient function as the correspondant value of G at the index multiplied by the defined gradient function
                    Gradienty_function_index = [Gradienty_function_index,m]; % Add the 'm' index in the temporary gradient function index
                    Physics2D.Gradienty_compute_Index{m} = 1; % Setting the computation gradient index
                % ELSEIF the value of the matrix at the index is zero or
                % the gradient function is not defined
                else
                    Physics2D.Gradienty_function{n,m} = @(X,Y) 0; % Storing zero in the gradient function matrix
                end
        end
        Physics2D.Gradienty_function_Index{n} = Gradienty_function_index; % Store the gradient function index for the 'm' index
    end
% ELSEIF there are 3 inputs and the gradient function is not defined
elseif (nargin == 3) && (iscell(Gradienty) == 0) && (isempty(Gradienty) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Gradienty_function_index = []; % Initializing the temporary gradient function index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF it is a diagonal term
                if (n == m)
                    Physics2D.Gradienty_function{n,m} = @(X,Y)Default_Gradienty(X,Y); % Storing the gradient function as the default gradient function
                    Gradienty_function_index = [Gradienty_function_index,m]; % Add the 'm' index in the temporary gradient function index
                    Physics2D.Gradienty_compute_Index{m} = 1; % Setting the computation gradient index
                % ELSEIF it is an extradiagonal term
                elseif (n ~= m)
                    Physics2D.Gradienty_function{n,m} = @(X,Y)0; % Storing zero in the gradient function matrix
                end
        end
        Physics2D.Gradienty_function_Index{n} = Gradienty_function_index; % Store the gradient function index for the 'm' index
    end
% ELSEIF there are 3 inputs and the gradient function is defined but not a
% cell array
elseif (nargin == 3) && (iscell(Gradienty) == 0) && (isempty(Gradienty) == 0)
    % FOR each component
    for n = 1:Method.Ncomponents
        Gradienty_function_index = []; % Initializing the temporary gradient function index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF it is a diagonal term
                if (n == m)
                    Physics2D.Gradienty_function{n,m} = @(X,Y)Gradienty(X,Y); % Storing the gradient function as the defined gradient function
                    Gradienty_function_index = [Gradienty_function_index,m]; % Add the 'm' index in the temporary gradient function index
                    Physics2D.Gradienty_compute_Index{m} = 1; % Setting the computation gradient index
                % ELSEIF it is an extradiagonal term
                elseif (n ~= m)
                    Physics2D.Gradienty_function{n,m} = @(X,Y)0; % Storing zero in the gradient function matrix
                end
        end
        Physics2D.Gradienty_function_Index{n} = Gradienty_function_index; % Store the gradient function index for the 'm' index
    end
% ELSEIF there are 3 inputs and the gradient function is a cell array
elseif (nargin == 3) && (iscell(Gradienty) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Gradienty_function_index = []; % Initializing the temporary gradient function index
        % FOR each component
        for m = 1:Method.Ncomponents
            % IF the gradient function is defined
            if (isempty(Gradienty{n,m}) == 0)
                Physics2D.Gradienty_function{n,m} = @(X,Y)Gradienty{n,m}(X,Y); % Storing the gradient function as the defined gradient function
                Gradienty_function_index = [Gradienty_function_index,m]; % Add the 'm' index in the temporary potential index
                Physics2D.Gradienty_compute_Index{m} = 1; % Setting the computation gradient index
            %ELSE if the gradient function is not defined
            else
                Physics2D.Gradienty_function{n,m} = @(X,Y)0; % Storing zero in the gradient function matrix
            end
        end
        Physics2D.Gradienty_function_Index{n} = Gradienty_function_index; % Store the gradient function index for the 'm' index
    end
% ELSEIF there are 2 inputs
elseif (nargin == 2)
    % FOR each component
    for n = 1:Method.Ncomponents
        Gradienty_function_index = []; % Initializing the temporary gradient function index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF it is a diagonal term
                if (n == m)
                    Physics2D.Gradienty_function{n,m} = @(X,Y)Default_Gradienty(X,Y); % Storing the gradient function as the default gradient function
                    Gradienty_function_index = [Gradienty_function_index,m]; % Add the 'm' index in the gradient function potential index
                    Physics2D.Gradienty_compute_Index{m} = 1; % Setting the computation gradient index
                % ELSEIF it is an extradiagonal term
                elseif (n ~= m)
                    Physics2D.Gradienty_function{n,m} = @(X,Y)0; % Storing zero in the gradient function matrix
                end
        end
        Physics2D.Gradienty_function_Index{n} = Gradienty_function_index; % Store the gradient function index for the 'm' index
    end
end