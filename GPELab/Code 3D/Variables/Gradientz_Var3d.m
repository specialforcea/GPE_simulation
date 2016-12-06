%% Addition of the gradient in the z direction function matrix in the physics of the
%% problem
%% INPUT:
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var3d.m)
%%          Physics3D: Structure containing variables concerning the physics of the problem in 3D (structure) (see Physics3D_Var3d.m)
%% INPUTS(OPTIONAL):
%%          Gradientz: Function or cell array of functions that depends on (X,Y,Z) which will be multiplied by the gradient in the z direction of the wave function in the physical problem (function or cell array of functions)
%%          (In the case of a function, the function will be applied only on the diagonal terms)
%%          G: Matrix that will be multiplied element by element to the gradient in the z direction function matrix (matrix)
%% OUTPUT:
%%          Physics3D: Structure containing variables concerning the physics of the problem in 3D (structure) (see Physics3D_Var3d.m)

function [Physics3D] = Gradientz_Var3d(Method, Physics3D, Gradientz, G)
%% Initializing the default gradient function
Default_Gradientz = @(X,Y,Z)(-1i)*(Physics3D.Omega(1)*Y-Physics3D.Omega(2)*X);

%% Adding the gradient in the y direction function matrix
% IF there are 4 inputs and the gradient function is not defined
if (nargin == 4) && (iscell(Gradientz) == 0) && (isempty(Gradientz) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Gradientz_function_index = []; % Initializing the temporary gradient index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero
                if (G(n,m)~= 0)
                    % IF it is a diagonal term
                    if(n == m)
                    Physics3D.Gradientz_function{n,m} = @(X,Y,Z)G(n,m)*Default_Gradientz(X,Y,Z); % Storing the gradient function as the correspondant value of G at the index added by the default gradient function
                    Gradientz_function_index = [Gradientz_function_index,m]; % Add the 'm' index in the temporary gradient function index
                    Physics3D.Gradientz_compute_Index{m} = 1; % Setting the computation gradient index
                    % ELSE if it is an extradiagonal term
                    else
                    Physics3D.Gradientz_function{n,m} = @(X,Y,Z) 0; % Storing zero in the gradient function matrix
                    end
                % ELSEIF the value of the matrix at the index is zero
                elseif (G(n,m) == 0)
                    Physics3D.Gradientz_function{n,m} = @(X,Y,Z) 0; % Storing zero in the gradient function matrix
                end
        end
        Physics3D.Gradientz_function_Index{n} = Gradientz_function_index; % Store the gradient function index for the 'm' index
    end
% ELSEIF there are 4 inputs and the gradient function is defined but not a
% cell array
elseif (nargin == 4) && (iscell(Gradientz) == 0) && (isempty(Gradientz) == 0)
    % FOR each component
    for n = 1:Method.Ncomponents
        Gradientz_function_index = []; % Initializing the temporary gradient function index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero
                if (G(n,m)~= 0)
                    Physics3D.Gradientz_function{n,m} = @(X,Y,Z)G(n,m)*Gradientz(X,Y,Z); % Storing the gradient function as the correspondant value of G at the index multiplied by the defined gradient function
                    Gradientz_function_index = [Gradientz_function_index,m]; % Add the 'm' index in the temporary gradient function index
                    Physics3D.Gradientz_compute_Index{m} = 1; % Setting the computation gradient index
                % ELSEIF the value of the matrix at the index is zero
                elseif (G(n,m) == 0)
                    Physics3D.Gradientz_function{n,m} = @(X,Y,Z)0; % Storing zero in the gradient function matrix
                end
        end
        Physics3D.Gradientz_function_Index{n} = Gradientz_function_index; % Store the gradient function index for the 'm' index
    end
% ELSEIF there are 4 inputs and the gradient function is a cell array
elseif (nargin == 4) && (iscell(Gradientz) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Gradientz_function_index = []; % Initializing the temporary gradient function index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero and
                % the gradient function is defined
                if (G(n,m)~= 0) && (isempty(Gradientz{n,m}) == 0)
                    Physics3D.Gradientz_function{n,m} = @(X,Y,Z)G(n,m)*Gradientz{n,m}(X,Y,Z); % Storing the gradient function as the correspondant value of G at the index multiplied by the defined gradient function
                    Gradientz_function_index = [Gradientz_function_index,m]; % Add the 'm' index in the temporary gradient function index
                    Physics3D.Gradientz_compute_Index{m} = 1; % Setting the computation gradient index
                % ELSEIF the value of the matrix at the index is zero or
                % the gradient function is not defined
                else
                    Physics3D.Gradientz_function{n,m} = @(X,Y,Z) 0; % Storing zero in the gradient function matrix
                end
        end
        Physics3D.Gradientz_function_Index{n} = Gradientz_function_index; % Store the gradient function index for the 'm' index
    end
% ELSEIF there are 3 inputs and the gradient function is not defined
elseif (nargin == 3) && (iscell(Gradientz) == 0) && (isempty(Gradientz) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Gradientz_function_index = []; % Initializing the temporary gradient function index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF it is a diagonal term
                if (n == m)
                    Physics3D.Gradientz_function{n,m} = @(X,Y,Z)Default_Gradientz(X,Y,Z); % Storing the gradient function as the default gradient function
                    Gradientz_function_index = [Gradientz_function_index,m]; % Add the 'm' index in the temporary gradient function index
                    Physics3D.Gradientz_compute_Index{m} = 1; % Setting the computation gradient index
                % ELSEIF it is an extradiagonal term
                elseif (n ~= m)
                    Physics3D.Gradientz_function{n,m} = @(X,Y,Z)0; % Storing zero in the gradient function matrix
                end
        end
        Physics3D.Gradientz_function_Index{n} = Gradientz_function_index; % Store the gradient function index for the 'm' index
    end
% ELSEIF there are 3 inputs and the gradient function is defined but not a
% cell array
elseif (nargin == 3) && (iscell(Gradientz) == 0) && (isempty(Gradientz) == 0)
    % FOR each component
    for n = 1:Method.Ncomponents
        Gradientz_function_index = []; % Initializing the temporary gradient function index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF it is a diagonal term
                if (n == m)
                    Physics3D.Gradientz_function{n,m} = @(X,Y,Z)Gradientz(X,Y,Z); % Storing the gradient function as the defined gradient function
                    Gradientz_function_index = [Gradientz_function_index,m]; % Add the 'm' index in the temporary gradient function index
                    Physics3D.Gradientz_compute_Index{m} = 1; % Setting the computation gradient index
                % ELSEIF it is an extradiagonal term
                elseif (n ~= m)
                    Physics3D.Gradientz_function{n,m} = @(X,Y,Z)0; % Storing zero in the gradient function matrix
                end
        end
        Physics3D.Gradientz_function_Index{n} = Gradientz_function_index; % Store the gradient function index for the 'm' index
    end
% ELSEIF there are 3 inputs and the gradient function is a cell array
elseif (nargin == 3) && (iscell(Gradientz) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Gradientz_function_index = []; % Initializing the temporary gradient function index
        % FOR each component
        for m = 1:Method.Ncomponents
            % IF the gradient function is defined
            if (isempty(Gradientz{n,m}) == 0)
                Physics3D.Gradientz_function{n,m} = @(X,Y,Z)Gradientz{n,m}(X,Y,Z); % Storing the gradient function as the defined gradient function
                Gradientz_function_index = [Gradientz_function_index,m]; % Add the 'm' index in the temporary potential index
                Physics3D.Gradientz_compute_Index{m} = 1; % Setting the computation gradient index
            %ELSE if the gradient function is not defined
            else
                Physics3D.Gradientz_function{n,m} = @(X,Y,Z)0; % Storing zero in the gradient function matrix
            end
        end
        Physics3D.Gradientz_function_Index{n} = Gradientz_function_index; % Store the gradient function index for the 'm' index
    end
% ELSEIF there are 2 inputs
elseif (nargin == 2)
    % FOR each component
    for n = 1:Method.Ncomponents
        Gradientz_function_index = []; % Initializing the temporary gradient function index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF it is a diagonal term
                if (n == m)
                    Physics3D.Gradientz_function{n,m} = @(X,Y,Z)Default_Gradientz(X,Y,Z); % Storing the gradient function as the default gradient function
                    Gradientz_function_index = [Gradientz_function_index,m]; % Add the 'm' index in the gradient function potential index
                    Physics3D.Gradientz_compute_Index{m} = 1; % Setting the computation gradient index
                % ELSEIF it is an extradiagonal term
                elseif (n ~= m)
                    Physics3D.Gradientz_function{n,m} = @(X,Y,Z)0; % Storing zero in the gradient function matrix
                end
        end
        Physics3D.Gradientz_function_Index{n} = Gradientz_function_index; % Store the gradient function index for the 'm' index
    end
end