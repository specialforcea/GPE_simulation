%% Addition of the gradient in the x direction function matrix in the physics of the
%% problem
%% INPUTS:
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var3d.m)
%%          Physics3D: Structure containing variables concerning the physics of the problem in 3D (structure) (see Physics3D_Var3d.m)
%% INPUTS(OPTIONAL):
%%          Gradientx: Function or cell array of functions that depends on (X,Y,Z) which will be multiplied by the gradient in the x direction of the wave function in the physical problem (function or cell array of functions)
%%          (In the case of a function, the function will be applied only on the diagonal terms)
%%          G: Matrix that will be multiplied element by element to the gradient in the x direction function matrix (matrix)
%% OUTPUT:
%%          Physics3D: Structure containing variables concerning the physics of the problem in 3D (structure) (see Physics3D_Var3d.m)

function [Physics3D] = Gradientx_Var3d(Method, Physics3D, Gradientx, G)
%% Initializing the default gradient function
Default_Gradientx = @(X,Y,Z)(-1i)*(Physics3D.Omega(2)*Z - Physics3D.Omega(3)*Y);

%% Adding the gradient in the y direction function matrix
% IF there are 4 inputs and the gradient function is not defined
if (nargin == 4) && (iscell(Gradientx) == 0) && (isempty(Gradientx) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Gradientx_function_index = []; % Initializing the temporary gradient index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero
                if (G(n,m)~= 0)
                    % IF it is a diagonal term
                    if(n == m)
                    Physics3D.Gradientx_function{n,m} = @(X,Y,Z)G(n,m)*Default_Gradientx(X,Y,Z); % Storing the gradient function as the correspondant value of G at the index added by the default gradient function
                    Gradientx_function_index = [Gradientx_function_index,m]; % Add the 'm' index in the temporary gradient function index
                    Physics3D.Gradientx_compute_Index{m} = 1; % Setting the computation gradient index
                    % ELSE if it is an extradiagonal term
                    else
                    Physics3D.Gradientx_function{n,m} = @(X,Y,Z) 0; % Storing zero in the gradient function matrix
                    end
                % ELSEIF the value of the matrix at the index is zero
                elseif (G(n,m) == 0)
                    Physics3D.Gradientx_function{n,m} = @(X,Y,Z) 0; % Storing zero in the gradient function matrix
                end
        end
        Physics3D.Gradientx_function_Index{n} = Gradientx_function_index; % Store the gradient function index for the 'm' index
    end
% ELSEIF there are 4 inputs and the gradient function is defined but not a
% cell array
elseif (nargin == 4) && (iscell(Gradientx) == 0) && (isempty(Gradientx) == 0)
    % FOR each component
    for n = 1:Method.Ncomponents
        Gradientx_function_index = []; % Initializing the temporary gradient function index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero
                if (G(n,m)~= 0)
                    Physics3D.Gradientx_function{n,m} = @(X,Y,Z)G(n,m)*Gradientx(X,Y,Z); % Storing the gradient function as the correspondant value of G at the index multiplied by the defined gradient function
                    Gradientx_function_index = [Gradientx_function_index,m]; % Add the 'm' index in the temporary gradient function index
                    Physics3D.Gradientx_compute_Index{m} = 1; % Setting the computation gradient index
                % ELSE if the value of the matrix at the index is zero
                else
                    Physics3D.Gradientx_function{n,m} = @(X,Y,Z) 0; % Storing zero in the gradient function matrix
                end
        end
        Physics3D.Gradientx_function_Index{n} = Gradientx_function_index; % Store the gradient function index for the 'm' index
    end
% ELSEIF there are 4 inputs and the gradient function is a cell array
elseif (nargin == 4) && (iscell(Gradientx) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Gradientx_function_index = []; % Initializing the temporary gradient function index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero and
                % the gradient function is defined
                if (G(n,m)~= 0) && (isempty(Gradientx{n,m}) == 0)
                    Physics3D.Gradientx_function{n,m} = @(X,Y,Z)G(n,m)*Gradientx{n,m}(X,Y,Z); % Storing the gradient function as the correspondant value of G at the index multiplied by the defined gradient function
                    Gradientx_function_index = [Gradientx_function_index,m]; % Add the 'm' index in the temporary gradient function index
                    Physics3D.Gradientx_compute_Index{m} = 1; % Setting the computation gradient index
                % ELSEIF the value of the matrix at the index is zero or
                % the gradient function is not defined
                else
                    Physics3D.Gradientx_function{n,m} = @(X,Y,Z) 0; % Storing zero in the gradient function matrix
                end
        end
        Physics3D.Gradientx_function_Index{n} = Gradientx_function_index; % Store the gradient function index for the 'm' index
    end
% ELSEIF there are 3 inputs and the gradient function is not defined
elseif (nargin == 3) && (iscell(Gradientx) == 0) && (isempty(Gradientx) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Gradientx_function_index = []; % Initializing the temporary gradient function index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF it is a diagonal term
                if (n == m)
                    Physics3D.Gradientx_function{n,m} = @(X,Y,Z)Default_Gradientx(X,Y,Z); % Storing the gradient function as the default gradient function
                    Gradientx_function_index = [Gradientx_function_index,m]; % Add the 'm' index in the temporary gradient function index
                    Physics3D.Gradientx_compute_Index{m} = 1; % Setting the computation gradient index
                % ELSEIF it is an extradiagonal term
                elseif (n ~= m)
                    Physics3D.Gradientx_function{n,m} = @(X,Y,Z)0; % Storing zero in the gradient function matrix
                end
        end
        Physics3D.Gradientx_function_Index{n} = Gradientx_function_index; % Store the gradient function index for the 'm' index
    end
% ELSEIF there are 3 inputs and the gradient function is defined but not a
% cell array
elseif (nargin == 3) && (iscell(Gradientx) == 0) && (isempty(Gradientx) == 0)
    % FOR each component
    for n = 1:Method.Ncomponents
        Gradientx_function_index = []; % Initializing the temporary gradient function index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF it is a diagonal term
                if (n == m)
                    Physics3D.Gradientx_function{n,m} = @(X,Y,Z)Gradientx(X,Y,Z); % Storing the gradient function as the defined gradient function
                    Gradientx_function_index = [Gradientx_function_index,m]; % Add the 'm' index in the temporary gradient function index
                    Physics3D.Gradientx_compute_Index{m} = 1; % Setting the computation gradient index
                % ELSEIF it is an extradiagonal term
                elseif (n ~= m)
                    Physics3D.Gradientx_function{n,m} = @(X,Y,Z)0; % Storing zero in the gradient function matrix
                end
        end
        Physics3D.Gradientx_function_Index{n} = Gradientx_function_index; % Store the gradient function index for the 'm' index
    end
% ELSEIF there are 3 inputs and the gradient function is a cell array
elseif (nargin == 3) && (iscell(Gradientx) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Gradientx_function_index = []; % Initializing the temporary gradient function index
        % FOR each component
        for m = 1:Method.Ncomponents
            % IF the gradient function is defined
            if (isempty(Gradientx{n,m}) == 0)
                Physics3D.Gradientx_function{n,m} = @(X,Y,Z)Gradientx{n,m}(X,Y,Z); % Storing the gradient function as the defined gradient function
                Gradientx_function_index = [Gradientx_function_index,m]; % Add the 'm' index in the temporary potential index
                Physics3D.Gradientx_compute_Index{m} = 1; % Setting the computation gradient index
            %ELSE if the gradient function is not defined
            else
                Physics3D.Gradientx_function{n,m} = @(X,Y,Z)0; % Storing zero in the gradient function matrix
            end
        end
        Physics3D.Gradientx_function_Index{n} = Gradientx_function_index; % Store the gradient function index for the 'm' index
    end
% ELSEIF there are 2 inputs
elseif (nargin == 2)
    % FOR each component
    for n = 1:Method.Ncomponents
        Gradientx_function_index = []; % Initializing the temporary gradient function index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF it is a diagonal term
                if (n == m)
                    Physics3D.Gradientx_function{n,m} = @(X,Y,Z)Default_Gradientx(X,Y,Z); % Storing the gradient function as the default gradient function
                    Gradientx_function_index = [Gradientx_function_index,m]; % Add the 'm' index in the gradient function potential index
                    Physics3D.Gradientx_compute_Index{m} = 1; % Setting the computation gradient index
                % ELSEIF it is an extradiagonal term
                elseif (n ~= m)
                    Physics3D.Gradientx_function{n,m} = @(X,Y,Z)0; % Storing zero in the gradient function matrix
                end
        end
        Physics3D.Gradientx_function_Index{n} = Gradientx_function_index; % Store the gradient function index for the 'm' index
    end
end