%% Addition of the dispersion function matrix in the physics of the
%% problem
%% INPUTS:
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var1d.m)
%%          Physics1D: Structure containing variables concerning the physics of the problem in 1D (structure) (see Physics1D_Var1d.m)
%% INPUTS(OPTIONAL):
%%          Dispersion: Function or cell array of functions that depends on (t,FFTX) which will be multiplied by the wave function in the physical problem (function or cell array of functions)
%%          (In the case of a function, the function will be applied only on the diagonal terms)
%%          G: Matrix that will be multiplied element by element to the potential function matrix (matrix)
%% OUTPUT:
%%          Physics1D: Structure containing variables concerning the physics of the problem in 1D (structure) (see Physics1D_Var1d.m)

function [Physics1D] = TimeDispersion_Var1d(Method, Physics1D, TimeDispersion, G)
%% Initializing the default potential
Default_Dispersion = @(t,FFTX) Physics1D.Delta*(FFTX.^2);

%% Adding the potential function matrix
% IF there are 4 inputs and the potential is not defined
if (nargin == 4) && (iscell(TimeDispersion) == 0) && (isempty(TimeDispersion) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        TimeDispersion_function_index = []; % Initializing the temporary potential index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero
                if (G(n,m)~= 0)
                    % IF it is a diagonal term
                    if(n == m)
                    Physics1D.TimeDispersion_function{n,m} = @(t,FFTX)G(n,m)*Default_Dispersion(FFTX); % Storing the potential as the correspondant value of G at the index added by the default potential
                    TimeDispersion_function_index = [TimeDispersion_function_index,m]; % Add the 'm' index in the temporary potential index
                    elseif (n ~= m)
                    Physics1D.TimeDispersion_function{n,m} = @(t,FFTX) 0; % Storing zero in the potential function matrix
                    end
                % ELSEIF the value of the matrix at the index is zero
                elseif (G(n,m) == 0)
                    Physics1D.TimeDispersion_function{n,m} = @(t,FFTX) 0; % Storing zero in the potential function matrix
                end
        end
        Physics1D.TimeDispersion_function_Index{n} = TimeDispersion_function_index; % Store the potential index for the 'm' index
    end
% ELSEIF there are 4 inputs and the potential is defined but not a
% cell array
elseif (nargin == 4) && (iscell(TimeDispersion) == 0) && (isempty(TimeDispersion) == 0)
    % FOR each component
    for n = 1:Method.Ncomponents
        TimeDispersion_function_index = []; % Initializing the temporary potential index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero
                if (G(n,m)~= 0)
                    Physics1D.TimeDispersion_function{n,m} = @(t,FFTX)G(n,m)*TimeDispersion(t,FFTX); % Storing the potential as the correspondant value of G at the index multiplied by the defined potential
                    TimeDispersion_function_index = [TimeDispersion_function_index,m]; % Add the 'm' index in the temporary potential index
                % ELSEIF the value of the matrix at the index is zero
                elseif (G(n,m) == 0)
                    Physics1D.TimeDispersion_function{n,m} = @(t,FFTX)0; % Storing zero in the potential function matrix
                end
        end
        Physics1D.TimeDispersion_function_Index{n} = TimeDispersion_function_index; % Store the potential index for the 'm' index
    end
% ELSEIF there are 4 inputs and the potential is a cell array
elseif (nargin == 4) && (iscell(TimeDispersion) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        TimeDispersion_function_index = []; % Initializing the temporary potential index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero and
                % the potential function is defined
                if (G(n,m)~= 0) && (isempty(TimeDispersion{n,m}) == 0)
                    Physics1D.TimeDispersion_function{n,m} = @(t,FFTX)G(n,m)*TimeDispersion{n,m}(t,FFTX); % Storing the potential as the correspondant value of G at the index multiplied by the defined potential
                    TimeDispersion_function_index = [TimeDispersion_function_index,m]; % Add the 'm' index in the temporary potential index
                % ELSEIF the value of the matrix at the index is zero or
                % the potential function is not defined
                else
                    Physics1D.TimeDispersion_function{n,m} = @(t,FFTX)0; % Storing zero in the potential function matrix
                end
        end
        Physics1D.Dispersion_function_Index{n} = TimeDispersion_function_index; % Store the potential index for the 'm' index
    end
% ELSEIF there are 3 inputs and the potential is not defined
elseif (nargin == 3) && (iscell(TimeDispersion) == 0) && (isempty(TimeDispersion) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        TimeDispersion_function_index = []; % Initializing the temporary potential index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF it is a diagonal term
                if (n == m)
                    Physics1D.TimeDispersion_function{n,m} = @(t,FFTX)Default_Dispersion(t,FFTX); % Storing the potential as the default potential
                    TimeDispersion_function_index = [TimeDispersion_function_index,m]; % Add the 'm' index in the temporary potential index
                % ELSEIF it is an extradiagonal term
                elseif (n ~= m)
                    Physics1D.TimeDispersion_function{n,m} = @(t,FFTX)0; % Storing zero in the potential function matrix
                end
        end
        Physics1D.TimeDispersion_function_Index{n} = TimeDispersion_function_index; % Store the potential index for the 'm' index
    end
% ELSEIF there are 3 inputs and the potential is defined but not a
% cell array
elseif (nargin == 3) && (iscell(TimeDispersion) == 0) && (isempty(TimeDispersion) == 0)
    % FOR each component
    for n = 1:Method.Ncomponents
        TimeDispersion_function_index = []; % Initializing the temporary potential index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF it is a diagonal term
                if (n == m)
                    Physics1D.TimeDispersion_function{n,m} = @(t,FFTX)TimeDispersion(t,FFTX); % Storing the potential as the defined potential
                    TimeDispersion_function_index = [TimeDispersion_function_index,m]; % Add the 'm' index in the temporary potential index
                % ELSEIF it is an extradiagonal term
                elseif (n ~= m)
                    Physics1D.TimeDispersion_function{n,m} = @(t,FFTX)0; % Storing zero in the potential function matrix
                end
        end
        Physics1D.TimeDispersion_function_Index{n} = TimeDispersion_function_index; % Store the potential index for the 'm' index
    end
% ELSEIF there are 3 inputs and the potential is a cell array
elseif (nargin == 3) && (iscell(TimeDispersion) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        TimeDispersion_function_index = []; % Initializing the temporary potential index
        % FOR each component
        for m = 1:Method.Ncomponents
            % IF the gradient function is defined
            if (isempty(TimeDispersion{n,m}) == 0)
                Physics1D.Dispersion_function{n,m} = @(FFTX)TimeDispersion{n,m}(FFTX); % Storing the potential as the defined potential
                TimeDispersion_function_index = [TimeDispersion_function_index,m]; % Add the 'm' index in the temporary potential index
            %ELSE if the gradient function is not defined
            else
                Physics1D.Dispersion_function{n,m} = @(FFTX) 0; % Storing zero in the potential function matrix
            end
        end
        Physics1D.Dispersion_function_Index{n} = TimeDispersion_function_index; % Store the potential index for the 'm' index
    end
% ELSEIF there are 2 inputs
elseif (nargin == 2)
    % FOR each component
    for n = 1:Method.Ncomponents
        TimeDispersion_function_index = []; % Initializing the temporary potential index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF it is a diagonal term
                if (n == m)
                    Physics1D.TimeDispersion_function{n,m} = @(t,FFTX)Default_Dispersion(t,FFTX); % Storing the potential as the default potential
                    TimeDispersion_function_index = [TimeDispersion_function_index,m]; % Add the 'm' index in the temporary potential index
                % ELSEIF it is an extradiagonal term
                elseif (n ~= m)
                    Physics1D.TimeDispersion_function{n,m} = @(t,FFTX)0; % Storing zero in the potential function matrix
                end
        end
        Physics1D.TimeDispersion_function_Index{n} = TimeDispersion_function_index; % Store the potential index for the 'm' index
    end
end