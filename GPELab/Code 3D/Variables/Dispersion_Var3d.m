%% Addition of the dispersion function matrix in the physics of the
%% problem
%% INPUTS:
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var3d.m)
%%          Physics3D: Structure containing variables concerning the physics of the problem in 3D (structure) (see Physics3D_Var3d.m)
%% INPUTS(OPTIONAL):
%%          Dispersion: Function or cell array of functions that depends on (FFTX,FFTY,FFTZ) which will be multiplied by the wave function in the physical problem (function or cell array of functions)
%%          (In the case of a function, the function will be applied only on the diagonal terms)
%%          G: Matrix that will be multiplied element by element to the potential function matrix (matrix)
%% OUTPUT:
%%          Physics3D: Structure containing variables concerning the physics of the problem in 3D (structure) (see Physics3D_Var3d.m)

function [Physics3D] = Dispersion_Var3d(Method, Physics3D, Dispersion, G)
%% Initializing the default potential
Default_Dispersion = @(FFTX,FFTY,FFTZ) Physics3D.Delta*(FFTX.^2+FFTY.^2+FFTZ.^2);

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
                    Physics3D.Dispersion_function{n,m} = @(FFTX,FFTY,FFTZ)G(n,m)*Default_Dispersion(FFTX,FFTY,FFTZ); % Storing the potential as the correspondant value of G at the index added by the default potential
                    Dispersion_function_index = [Dispersion_function_index,m]; % Add the 'm' index in the temporary potential index
                    elseif (n ~= m)
                    Physics3D.Dispersion_function{n,m} = @(FFTX,FFTY,FFTZ) 0; % Storing zero in the potential function matrix
                    end
                % ELSEIF the value of the matrix at the index is zero
                elseif (G(n,m) == 0)
                    Physics3D.Dispersion_function{n,m} = @(FFTX,FFTY,FFTZ) 0; % Storing zero in the potential function matrix
                end
        end
        Physics3D.Dispersion_function_Index{n} = Dispersion_function_index; % Store the potential index for the 'm' index
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
                    Physics3D.Dispersion_function{n,m} = @(FFTX,FFTY,FFTZ)G(n,m)*Dispersion(FFTX,FFTY,FFTZ); % Storing the potential as the correspondant value of G at the index multiplied by the defined potential
                    Dispersion_function_index = [Dispersion_function_index,m]; % Add the 'm' index in the temporary potential index
                % ELSEIF the value of the matrix at the index is zero
                elseif (G(n,m) == 0)
                    Physics3D.Dispersion_function{n,m} = @(FFTX,FFTY,FFTZ)0; % Storing zero in the potential function matrix
                end
        end
        Physics3D.Dispersion_function_Index{n} = Dispersion_function_index; % Store the potential index for the 'm' index
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
                    Physics3D.Dispersion_function{n,m} = @(FFTX,FFTY,FFTZ)G(n,m)*Dispersion{n,m}(FFTX,FFTY,FFTZ); % Storing the potential as the correspondant value of G at the index multiplied by the defined potential
                    Dispersion_function_index = [Dispersion_function_index,m]; % Add the 'm' index in the temporary potential index
                % ELSEIF the value of the matrix at the index is zero or
                % the potential function is not defined
                else
                    Physics3D.Dispersion_function{n,m} = @(FFTX,FFTY,FFTZ)0; % Storing zero in the potential function matrix
                end
        end
        Physics3D.Dispersion_function_Index{n} = Dispersion_function_index; % Store the potential index for the 'm' index
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
                    Physics3D.Dispersion_function{n,m} = @(FFTX,FFTY,FFTZ)Default_Dispersion(FFTX,FFTY,FFTZ); % Storing the potential as the default potential
                    Dispersion_function_index = [Dispersion_function_index,m]; % Add the 'm' index in the temporary potential index
                % ELSEIF it is an extradiagonal term
                elseif (n ~= m)
                    Physics3D.Dispersion_function{n,m} = @(FFTX,FFTY,FFTZ)0; % Storing zero in the potential function matrix
                end
        end
        Physics3D.Dispersion_function_Index{n} = Dispersion_function_index; % Store the potential index for the 'm' index
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
                    Physics3D.Dispersion_function{n,m} = @(FFTX,FFTY,FFTZ)Dispersion(FFTX,FFTY,FFTZ); % Storing the potential as the defined potential
                    Dispersion_function_index = [Dispersion_function_index,m]; % Add the 'm' index in the temporary potential index
                % ELSEIF it is an extradiagonal term
                elseif (n ~= m)
                    Physics3D.Dispersion_function{n,m} = @(FFTX,FFTY,FFTZ)0; % Storing zero in the potential function matrix
                end
        end
        Physics3D.Dispersion_function_Index{n} = Dispersion_function_index; % Store the potential index for the 'm' index
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
                Physics3D.Dispersion_function{n,m} = @(FFTX,FFTY,FFTZ)Dispersion{n,m}(FFTX,FFTY,FFTZ); % Storing the potential as the defined potential
                Dispersion_function_index = [Dispersion_function_index,m]; % Add the 'm' index in the temporary potential index
            %ELSE if the gradient function is not defined
            else
                Physics3D.Dispersion_function{n,m} = @(FFTX,FFTY,FFTZ) 0; % Storing zero in the potential function matrix
            end
        end
        Physics3D.Dispersion_function_Index{n} = Dispersion_function_index; % Store the potential index for the 'm' index
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
                    Physics3D.Dispersion_function{n,m} = @(FFTX,FFTY,FFTZ)Default_Dispersion(FFTX,FFTY,FFTZ); % Storing the potential as the default potential
                    Dispersion_function_index = [Dispersion_function_index,m]; % Add the 'm' index in the temporary potential index
                % ELSEIF it is an extradiagonal term
                elseif (n ~= m)
                    Physics3D.Dispersion_function{n,m} = @(FFTX,FFTY,FFTZ)0; % Storing zero in the potential function matrix
                end
        end
        Physics3D.Dispersion_function_Index{n} = Dispersion_function_index; % Store the potential index for the 'm' index
    end
end