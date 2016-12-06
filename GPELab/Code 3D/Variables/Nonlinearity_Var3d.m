%% Addition of the nonlinearity function matrix in the physics of the
%% problem
%% INPUTS:
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var3d.m)
%%          Physics3D: Structure containing variables concerning the physics of the problem in 3D (structure) (see Physics3D_Var3d.m)
%% INPUTS(OPTIONAL):
%%          Nonlinearity: Function or cell array of functions that depends on (phi,X,Y,Z) which will be multiplied by the wave function in the physical problem (function or cell array of functions)
%%          (In the case of a function, the function will be applied only on the diagonal terms)
%%          G: Matrix that will be multiplied element by element to the nonlinearity function matrix (matrix)
%%          Nonlinearity_energy: Function or cell array of functions that depends on (phi,X,Y,Z) which will be multiplied by the wave function in the computation of the energy (function or cell array of functions)
%%          (In the case of a function, the function will be applied only on the diagonal terms. Moreover, if it is not defined, the energy will be computed using the nonlinear function used in the computation of the ground state)
%% OUTPUT:
%%          Physics3D: Structure containing variables concerning the physics of the problem in 3D (structure) (see Physics3D_Var3d.m)

function [Physics3D] = Nonlinearity_Var3d(Method, Physics3D, Nonlinearity, G, Nonlinearity_energy)
%% Initializing the nonlinearity index
Physics3D.Nonlinearity_function_Index = [];
%% Initializing the default nonlinearity
Default_Nonlinearity = Cubic3d(Method);
Default_Nonlinearity_energy = Cubic_energy3d(Method);

%% Adding the nonlinearity function matrix
% IF there are 5 inputs and the nonlinearity is not defined and G is not
% defined
if (nargin == 5) && (iscell(Nonlinearity) == 0) && (isempty(Nonlinearity) == 1) && (isempty(G) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Nonlinearity_function_index = []; % Initializing the temporary nonlinearity index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the nonlinearity at the index is defined
                if (isempty(Default_Nonlinearity{n,m}) == 0)
                    Physics3D.Nonlinearity_function{n,m} = @(Phi,X,Y,Z)Default_Nonlinearity{n,m}(Phi,X,Y,Z); % Storing the nonlinearity as the correspondant value of G at the index multiplied by the default nonlinearity
                    Physics3D.Nonlinearity_energy_function{n,m} = @(Phi,X,Y,Z)Nonlinearity_energy{n,m}(Phi,X,Y,Z); % Storing the energy's nonlinearity as the correspondant value of G at the index multiplied by the defined nonlinearity
                    Nonlinearity_function_index = [Nonlinearity_function_index,m]; % Add the 'm' index in the temporary nonlinearity index
                % ELSE if the value of the matrix at the index is zero or the
                % nonlinearity at the index is not defined
                else
                    Physics3D.Nonlinearity_function{n,m} = @(Phi,X,Y,Z)0; % Storing zero in the nonlinearity function matrix
                    Physics3D.Nonlinearity_energy_function{n,m} = @(Phi,X,Y,Z)0; % Storing zero in the energy's nonlinearity function matrix
                end
        end
        Physics3D.Nonlinearity_function_Index{n} = Nonlinearity_function_index; % Store the potential index for the 'm' index
    end
% ELSEIF there are 5 inputs and the nonlinearity is defined but not a
% cell array and G is not defined
elseif (nargin == 5) && (iscell(Nonlinearity) == 0) && (isempty(Nonlinearity) == 0) && (isempty(G) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Nonlinearity_function_index = []; % Initializing the temporary nonlinearity index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF it is a diagonal term
                if (m == n)
                    Physics3D.Nonlinearity_function{n,m} = @(Phi,X,Y,Z)Nonlinearity(Phi{n},X,Y,Z); % Storing the nonlinearity as the correspondant value of G at the index multiplied by the defined nonlinearity
                    Physics3D.Nonlinearity_energy_function{n,m} = @(Phi,X,Y,Z)Nonlinearity_energy(Phi{n},X,Y,Z); % Storing the energy's nonlinearity as the correspondant value of G at the index multiplied by the defined nonlinearity
                    Nonlinearity_function_index = [Nonlinearity_function_index,m]; % Add the 'm' index in the temporary nonlinearity index
                % ELSE if it is an extradiagonal term
                else
                    Physics3D.Nonlinearity_function{n,m} = @(Phi,X,Y,Z)0; % Storing zero in the nonlinearity function matrix
                    Physics3D.Nonlinearity_energy_function{n,m} = @(Phi,X,Y,Z)0; % Storing zero in the energy's nonlinearity function matrix
                end
        end
        Physics3D.Nonlinearity_function_Index{n} = Nonlinearity_function_index; % Store the nonlinearity index for the 'm' index
    end
% ELSEIF there are 5 inputs and the nonlinearity is a cell array and G is
% not defined
elseif (nargin == 5) && (iscell(Nonlinearity) == 1) && (isempty(G) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Nonlinearity_function_index = []; % Initializing the temporary nonlinearity index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero and
                % the nonlinearity at the index is defined
                if (isempty(Nonlinearity{n,m}) == 0)
                    Physics3D.Nonlinearity_function{n,m} = @(Phi,X,Y,Z)Nonlinearity{n,m}(Phi,X,Y,Z); % Storing the nonlinearity as the correspondant value of G at the index multiplied by the defined nonlinearity
                    Physics3D.Nonlinearity_energy_function{n,m} = @(Phi,X,Y,Z)Nonlinearity_energy{n,m}(Phi,X,Y,Z); % Storing the energy's nonlinearity as the correspondant value of G at the index multiplied by the defined nonlinearity
                    Nonlinearity_function_index = [Nonlinearity_function_index,m]; % Add the 'm' index in the temporary nonlinearity index
                % ELSE if the value of the matrix at the index is zero or the
                % nonlinearity at the index is not defined
                else
                    Physics3D.Nonlinearity_function{n,m} = @(Phi,X,Y,Z)0; % Storing zero in the nonlinearity function matrix
                    Physics3D.Nonlinearity_energy_function{n,m} = @(Phi,X,Y,Z)0; % Storing zero in the energy's nonlinearity function matrix
                end
        end
        Physics3D.Nonlinearity_function_Index{n} = Nonlinearity_function_index; % Store the nonlinearity index for the 'm' index
    end
% IF there are 5 inputs and the nonlinearity is not defined and G is
% defined
elseif (nargin == 5) && (iscell(Nonlinearity) == 0) && (isempty(Nonlinearity) == 1) && (isempty(G) == 0)
    % FOR each component
    for n = 1:Method.Ncomponents
        Nonlinearity_function_index = []; % Initializing the temporary nonlinearity index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero and
                % the nonlinearity at the index is defined
                if (G(n,m)~= 0) && (isempty(Default_Nonlinearity{n,m}) == 0)
                    Physics3D.Nonlinearity_function{n,m} = @(Phi,X,Y,Z)G(n,m)*Default_Nonlinearity{n,m}(Phi,X,Y,Z); % Storing the nonlinearity as the correspondant value of G at the index multiplied by the default nonlinearity
                    Physics3D.Nonlinearity_energy_function{n,m} = @(Phi,X,Y,Z)G(n,m)*Nonlinearity_energy{n,m}(Phi,X,Y,Z); % Storing the energy's nonlinearity as the correspondant value of G at the index multiplied by the defined nonlinearity
                    Nonlinearity_function_index = [Nonlinearity_function_index,m]; % Add the 'm' index in the temporary nonlinearity index
                % ELSE if the value of the matrix at the index is zero or the
                % nonlinearity at the index is not defined
                else
                    Physics3D.Nonlinearity_function{n,m} = @(Phi,X,Y,Z)0; % Storing zero in the nonlinearity function matrix
                    Physics3D.Nonlinearity_energy_function{n,m} = @(Phi,X,Y,Z)0; % Storing zero in the energy's nonlinearity function matrix
                end
        end
        Physics3D.Nonlinearity_function_Index{n} = Nonlinearity_function_index; % Store the potential index for the 'm' index
    end
% ELSEIF there are 5 inputs and the nonlinearity is defined but not a
% cell array and G is defined
elseif (nargin == 5) && (iscell(Nonlinearity) == 0) && (isempty(Nonlinearity) == 0) && (isempty(G) == 0)
    % FOR each component
    for n = 1:Method.Ncomponents
        Nonlinearity_function_index = []; % Initializing the temporary nonlinearity index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero
                if (G(n,m)~= 0)
                    Physics3D.Nonlinearity_function{n,m} = @(Phi,X,Y,Z)G(n,m)*Nonlinearity(Phi{n},X,Y,Z); % Storing the nonlinearity as the correspondant value of G at the index multiplied by the defined nonlinearity
                    Physics3D.Nonlinearity_energy_function{n,m} = @(Phi,X,Y,Z)G(n,m)*Nonlinearity_energy(Phi{n},X,Y,Z); % Storing the energy's nonlinearity as the correspondant value of G at the index multiplied by the defined nonlinearity
                    Nonlinearity_function_index = [Nonlinearity_function_index,m]; % Add the 'm' index in the temporary nonlinearity index
                % ELSEIF the value of the matrix at the index is zero
                elseif (G(n,m) == 0)
                    Physics3D.Nonlinearity_function{n,m} = @(Phi,X,Y,Z)0; % Storing zero in the nonlinearity function matrix
                    Physics3D.Nonlinearity_energy_function{n,m} = @(Phi,X,Y,Z)0; % Storing zero in the energy's nonlinearity function matrix
                end
        end
        Physics3D.Nonlinearity_function_Index{n} = Nonlinearity_function_index; % Store the nonlinearity index for the 'm' index
    end
% ELSEIF there are 5 inputs and the nonlinearity is a cell array and G is
% defined
elseif (nargin == 5) && (iscell(Nonlinearity) == 1) && (isempty(G) == 0)
    % FOR each component
    for n = 1:Method.Ncomponents
        Nonlinearity_function_index = []; % Initializing the temporary nonlinearity index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero and
                % the nonlinearity at the index is defined
                if (G(n,m)~= 0) && (isempty(Nonlinearity{n,m}) == 0)
                    Physics3D.Nonlinearity_function{n,m} = @(Phi,X,Y,Z)G(n,m)*Nonlinearity{n,m}(Phi,X,Y,Z); % Storing the nonlinearity as the correspondant value of G at the index multiplied by the defined nonlinearity
                    Physics3D.Nonlinearity_energy_function{n,m} = @(Phi,X,Y,Z)G(n,m)*Nonlinearity_energy(Phi{n},X,Y,Z); % Storing the energy's nonlinearity as the correspondant value of G at the index multiplied by the defined nonlinearity
                    Nonlinearity_function_index = [Nonlinearity_function_index,m]; % Add the 'm' index in the temporary nonlinearity index
                % ELSE if the value of the matrix at the index is zero or the
                % nonlinearity at the index is not defined
                else
                    Physics3D.Nonlinearity_function{n,m} = @(Phi,X,Y,Z)0; % Storing zero in the nonlinearity function matrix
                    Physics3D.Nonlinearity_energy_function{n,m} = @(Phi,X,Y,Z)0; % Storing zero in the energy's nonlinearity function matrix
                end
        end
        Physics3D.Nonlinearity_function_Index{n} = Nonlinearity_function_index; % Store the nonlinearity index for the 'm' index
    end
% IF there are 4 inputs and the nonlinearity is not defined
elseif (nargin == 4) && (iscell(Nonlinearity) == 0) && (isempty(Nonlinearity) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Nonlinearity_function_index = []; % Initializing the temporary nonlinearity index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero and
                % the nonlinearity at the index is defined
                if (G(n,m)~= 0) && (isempty(Default_Nonlinearity{n,m}) == 0)
                    Physics3D.Nonlinearity_function{n,m} = @(Phi,X,Y,Z)G(n,m)*Default_Nonlinearity{n,m}(Phi,X,Y,Z); % Storing the nonlinearity as the correspondant value of G at the index multiplied by the default nonlinearity
                    Physics3D.Nonlinearity_energy_function{n,m} = @(Phi,X,Y,Z)G(n,m)*Default_Nonlinearity_energy{n,m}(Phi,X,Y,Z); % Storing the energy's nonlinearity as the correspondant value of G at the index multiplied by the default nonlinearity
                    Nonlinearity_function_index = [Nonlinearity_function_index,m]; % Add the 'm' index in the temporary nonlinearity index
                % ELSE if the value of the matrix at the index is zero or the
                % nonlinearity at the index is not defined
                else
                    Physics3D.Nonlinearity_function{n,m} = @(Phi,X,Y,Z)0; % Storing zero in the nonlinearity function matrix
                end
        end
        Physics3D.Nonlinearity_function_Index{n} = Nonlinearity_function_index; % Store the potential index for the 'm' index
    end
% ELSEIF there are 4 inputs and the nonlinearity is defined but not a
% cell array
elseif (nargin == 4) && (iscell(Nonlinearity) == 0) && (isempty(Nonlinearity) == 0)
    % FOR each component
    for n = 1:Method.Ncomponents
        Nonlinearity_function_index = []; % Initializing the temporary nonlinearity index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero
                if (G(n,m)~= 0)
                    Physics3D.Nonlinearity_function{n,m} = @(Phi,X,Y,Z)G(n,m)*Nonlinearity(Phi{n},X,Y,Z); % Storing the nonlinearity as the correspondant value of G at the index multiplied by the defined nonlinearity
                    Nonlinearity_function_index = [Nonlinearity_function_index,m]; % Add the 'm' index in the temporary nonlinearity index
                % ELSE if it is an extradiagonal term or the value of the matrix at the index is zero
                else
                    Physics3D.Nonlinearity_function{n,m} = @(Phi,X,Y,Z)0; % Storing zero in the nonlinearity function matrix
                end
        end
        Physics3D.Nonlinearity_function_Index{n} = Nonlinearity_function_index; % Store the nonlinearity index for the 'm' index
    end
% ELSEIF there are 4 inputs and the nonlinearity is a cell array
elseif (nargin == 4) && (iscell(Nonlinearity) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Nonlinearity_function_index = []; % Initializing the temporary nonlinearity index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero and
                % the nonlinearity at the index is defined
                if (G(n,m)~= 0) && (isempty(Nonlinearity{n,m}) == 0)
                    Physics3D.Nonlinearity_function{n,m} = @(Phi,X,Y,Z)G(n,m)*Nonlinearity{n,m}(Phi,X,Y,Z); % Storing the nonlinearity as the correspondant value of G at the index multiplied by the defined nonlinearity
                    Nonlinearity_function_index = [Nonlinearity_function_index,m]; % Add the 'm' index in the temporary nonlinearity index
                % ELSE if the value of the matrix at the index is zero or the
                % nonlinearity at the index is not defined
                else
                    Physics3D.Nonlinearity_function{n,m} = @(Phi,X,Y,Z)0; % Storing zero in the nonlinearity function matrix
                end
        end
        Physics3D.Nonlinearity_function_Index{n} = Nonlinearity_function_index; % Store the nonlinearity index for the 'm' index
    end
% ELSEIF there are 3 inputs and the nonlinearity is not defined
elseif (nargin == 3) && (iscell(Nonlinearity) == 0) && (isempty(Nonlinearity) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Nonlinearity_function_index = []; % Initializing the temporary nonlinearity index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the nonlinearity is defined
                if (isempty(Default_Nonlinearity{n,m}) == 0)
                    Physics3D.Nonlinearity_function{n,m} = @(Phi,X,Y,Z)Default_Nonlinearity{n,m}(Phi,X,Y,Z); % Storing the nonlinearity as the default nonlinearity
                    Physics3D.Nonlinearity_energy_function{n,m} = @(Phi,X,Y,Z)Default_Nonlinearity_energy{n,m}(Phi,X,Y,Z); % Storing the energy's nonlinearity as the correspondant value of G at the index multiplied by the default nonlinearity
                    Nonlinearity_function_index = [Nonlinearity_function_index,m]; % Add the 'm' index in the temporary potential index
                % ELSE if the nonlinearity is not defined
                else
                    Physics3D.Nonlinearity_function{n,m} = @(Phi,X,Y,Z)0; % Storing zero in the nonlinearity function matrix
                end
        end
        Physics3D.Nonlinearity_function_Index{n} = Nonlinearity_function_index; % Store the nonlinearity index for the 'm' index
    end
% ELSEIF there are 3 inputs and the nonlinearity is defined but not a
% cell array
elseif (nargin == 3) && (iscell(Nonlinearity) == 0) && (isempty(Nonlinearity) == 0)
    % FOR each component
    for n = 1:Method.Ncomponents
        Nonlinearity_function_index = []; % Initializing the nonlinearity potential index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF it is a diagonal term
                if (n == m)
                    Physics3D.Nonlinearity_function{n,m} = @(Phi,X,Y,Z)Nonlinearity(Phi{n},X,Y,Z); % Storing the nonlinearity as the defined nonlinearity
                    Nonlinearity_function_index = [Nonlinearity_function_index,m]; % Add the 'm' index in the temporary nonlinearity index
                % ELSEIF it is an extradiagonal term
                elseif (n ~= m)
                    Physics3D.Nonlinearity_function{n,m} = @(Phi,X,Y,Z)0; % Storing zero in the nonlinearity function matrix
                end
        end
        Physics3D.Nonlinearity_function_Index{n} = Nonlinearity_function_index; % Store the nonlinearity index for the 'm' index
    end
% ELSEIF there are 3 inputs and the nonlinearity is a cell array
elseif (nargin == 3) && (iscell(Nonlinearity) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Nonlinearity_function_index = []; % Initializing the temporary nonlinearity index
        % FOR each component
        for m = 1:Method.Ncomponents
                    Physics3D.Nonlinearity_function{n,m} = @(Phi,X,Y,Z)Nonlinearity{n,m}(Phi,X,Y,Z); % Storing the nonlinearity as the defined nonlinearity
                    Nonlinearity_function_index = [Nonlinearity_function_index,m]; % Add the 'm' index in the temporary nonlinearity index
        end
        Physics3D.Nonlinearity_function_Index{n} = Nonlinearity_function_index; % Store the nonlinearity index for the 'm' index
    end
% ELSEIF there are 2 inputs
elseif (nargin == 2)
    % FOR each component
    for n = 1:Method.Ncomponents
        Nonlinearity_function_index = []; % Initializing the temporary nonlinearity index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF it is a diagonal term
                if (isempty(Default_Nonlinearity{n,m}) == 0)
                    Physics3D.Nonlinearity_function{n,m} = @(Phi,X,Y,Z)Default_Nonlinearity{n,m}(Phi,X,Y,Z); % Storing the nonlinearity as the default nonlinearity
                    Physics3D.Nonlinearity_energy_function{n,m} = @(Phi,X,Y,Z)Default_Nonlinearity_energy{n,m}(Phi,X,Y,Z); % Storing the energy's nonlinearity as the correspondant value of G at the index multiplied by the default nonlinearity
                    Nonlinearity_function_index = [Nonlinearity_function_index,m]; % Add the 'm' index in the temporary nonlinearity index
                % ELSEIF it is an extradiagonal term
                else
                    Physics3D.Nonlinearity{n,m} = @(Phi,X,Y,Z)0; % Storing zero in the nonlinearity function matrix
                end
        end
        Physics3D.Nonlinearity_function_Index{n} = Nonlinearity_function_index; % Store the nonlinearity index for the 'm' index
    end
end