%% Addition of the dispersion function matrix in the physics of the
%% problem
%% INPUTS:
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          Physics2D: Structure containing variables concerning the physics of the problem in 2D (structure) (see Physics2D_Var2d.m)
%% INPUTS(OPTIONAL):
%%          Dispersion: Function or cell array of functions that depends on (FFTX,FFTY) which will be multiplied by the wave function in the physical problem (function or cell array of functions)
%%          (In the case of a function, the function will be applied only on the diagonal terms)
%%          G: Matrix that will be multiplied element by element to the potential function matrix (matrix)
%% OUTPUT:
%%          Physics2D: Structure containing variables concerning the physics of the problem in 2D (structure) (see Physics2D_Var2d.m)

function [Physics2D] = Dispersion_Var2d(Method, Physics2D, Dispersion, G)
%% Initializing the default potential
Default_Dispersion = @(FFTX,FFTY) Physics2D.Delta*(FFTX.^2+FFTY.^2);

%% Adding the potential function matrix
% IF there are 4 inputs and the potential is not defined
if (nargin == 4) && (iscell(Dispersion) == 0) && (isempty(Dispersion) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Dispersion_function_index = []; % Initializing the temporary potential index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero
                if (G(n,m)~= 0)
                    % IF it is a diagonal term
                    if(n == m)
                    Physics2D.Dispersion_function{n,m} = @(FFTX,FFTY)G(n,m)*Default_Dispersion(FFTX,FFTY); % Storing the potential as the correspondant value of G at the index added by the default potential
                    Dispersion_function_index = [Dispersion_function_index,m]; % Add the 'm' index in the temporary potential index
                    elseif (n ~= m)
                    Physics2D.Dispersion_function{n,m} = @(FFTX,FFTY) 0; % Storing zero in the potential function matrix
                    end
                % ELSEIF the value of the matrix at the index is zero
                elseif (G(n,m) == 0)
                    Physics2D.Dispersion_function{n,m} = @(FFTX,FFTY) 0; % Storing zero in the potential function matrix
                end
        end
        Physics2D.Dispersion_function_Index{n} = Dispersion_function_index; % Store the potential index for the 'm' index
    end
% ELSEIF there are 4 inputs and the potential is defined but not a
% cell array
elseif (nargin == 4) && (iscell(Dispersion) == 0) && (isempty(Dispersion) == 0)
    % FOR each component
    for n = 1:Method.Ncomponents
        Dispersion_function_index = []; % Initializing the temporary potential index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero
                if (G(n,m)~= 0)
                    Physics2D.Dispersion_function{n,m} = @(FFTX,FFTY)G(n,m)*Dispersion(FFTX,FFTY); % Storing the potential as the correspondant value of G at the index multiplied by the defined potential
                    Dispersion_function_index = [Dispersion_function_index,m]; % Add the 'm' index in the temporary potential index
                % ELSEIF the value of the matrix at the index is zero
                elseif (G(n,m) == 0)
                    Physics2D.Dispersion_function{n,m} = @(FFTX,FFTY)0; % Storing zero in the potential function matrix
                end
        end
        Physics2D.Dispersion_function_Index{n} = Dispersion_function_index; % Store the potential index for the 'm' index
    end
% ELSEIF there are 4 inputs and the potential is a cell array
elseif (nargin == 4) && (iscell(Dispersion) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Dispersion_function_index = []; % Initializing the temporary potential index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero and
                % the potential function is defined
                if (G(n,m)~= 0) && (isempty(Dispersion{n,m}) == 0)
                    Physics2D.Dispersion_function{n,m} = @(FFTX,FFTY)G(n,m)*Dispersion{n,m}(FFTX,FFTY); % Storing the potential as the correspondant value of G at the index multiplied by the defined potential
                    Dispersion_function_index = [Dispersion_function_index,m]; % Add the 'm' index in the temporary potential index
                % ELSEIF the value of the matrix at the index is zero or
                % the potential function is not defined
                else
                    Physics2D.Dispersion_function{n,m} = @(FFTX,FFTY)0; % Storing zero in the potential function matrix
                end
        end
        Physics2D.Dispersion_function_Index{n} = Dispersion_function_index; % Store the potential index for the 'm' index
    end
% ELSEIF there are 3 inputs and the potential is not defined
elseif (nargin == 3) && (iscell(Dispersion) == 0) && (isempty(Dispersion) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Dispersion_function_index = []; % Initializing the temporary potential index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF it is a diagonal term
                if (n == m)
                    Physics2D.Dispersion_function{n,m} = @(FFTX,FFTY)Default_Dispersion(FFTX,FFTY); % Storing the potential as the default potential
                    Dispersion_function_index = [Dispersion_function_index,m]; % Add the 'm' index in the temporary potential index
                % ELSEIF it is an extradiagonal term
                elseif (n ~= m)
                    Physics2D.Dispersion_function{n,m} = @(FFTX,FFTY)0; % Storing zero in the potential function matrix
                end
        end
        Physics2D.Dispersion_function_Index{n} = Dispersion_function_index; % Store the potential index for the 'm' index
    end
% ELSEIF there are 3 inputs and the potential is defined but not a
% cell array
elseif (nargin == 3) && (iscell(Dispersion) == 0) && (isempty(Dispersion) == 0)
    % FOR each component
    for n = 1:Method.Ncomponents
        Dispersion_function_index = []; % Initializing the temporary potential index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF it is a diagonal term
                if (n == m)
                    Physics2D.Dispersion_function{n,m} = @(FFTX,FFTY)Dispersion(FFTX,FFTY); % Storing the potential as the defined potential
                    Dispersion_function_index = [Dispersion_function_index,m]; % Add the 'm' index in the temporary potential index
                % ELSEIF it is an extradiagonal term
                elseif (n ~= m)
                    Physics2D.Dispersion_function{n,m} = @(FFTX,FFTY)0; % Storing zero in the potential function matrix
                end
        end
        Physics2D.Dispersion_function_Index{n} = Dispersion_function_index; % Store the potential index for the 'm' index
    end
% ELSEIF there are 3 inputs and the potential is a cell array
elseif (nargin == 3) && (iscell(Dispersion) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Dispersion_function_index = []; % Initializing the temporary potential index
        % FOR each component
        for m = 1:Method.Ncomponents
            % IF the gradient function is defined
            if (isempty(Dispersion{n,m}) == 0)
                Physics2D.Dispersion_function{n,m} = @(FFTX,FFTY)Dispersion{n,m}(FFTX,FFTY); % Storing the potential as the defined potential
                Dispersion_function_index = [Dispersion_function_index,m]; % Add the 'm' index in the temporary potential index
            %ELSE if the gradient function is not defined
            else
                Physics2D.Dispersion_function{n,m} = @(FFTX,FFTY) 0; % Storing zero in the potential function matrix
            end
        end
        Physics2D.Dispersion_function_Index{n} = Dispersion_function_index; % Store the potential index for the 'm' index
    end
% ELSEIF there are 2 inputs
elseif (nargin == 2)
    % FOR each component
    for n = 1:Method.Ncomponents
        Dispersion_function_index = []; % Initializing the temporary potential index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF it is a diagonal term
                if (n == m)
                    Physics2D.Dispersion_function{n,m} = @(FFTX,FFTY)Default_Dispersion(FFTX,FFTY); % Storing the potential as the default potential
                    Dispersion_function_index = [Dispersion_function_index,m]; % Add the 'm' index in the temporary potential index
                % ELSEIF it is an extradiagonal term
                elseif (n ~= m)
                    Physics2D.Dispersion_function{n,m} = @(FFTX,FFTY)0; % Storing zero in the potential function matrix
                end
        end
        Physics2D.Dispersion_function_Index{n} = Dispersion_function_index; % Store the potential index for the 'm' index
    end
end