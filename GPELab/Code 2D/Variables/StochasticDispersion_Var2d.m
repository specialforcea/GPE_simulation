%% Addition of the stochastic dispersion in the physics of the
%% problem
%% INPUTS:
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          Physics2D: Structure containing variables concerning the physics of the problem in 2D (structure) (see Physics2D_Var2d.m)
%% INPUTS(OPTIONAL):
%%          Dispersion: Function or cell array of functions that depends on (Phi,W,FFTX) which will be applied to the wave function in the physical problem (function or cell array of functions)
%%          (In the case of a function, the function will be applied only on the diagonal terms)
%%          G: Matrix that will be multiplied element by element to the potential function matrix (matrix)
%%          StochasticProcess: Function that depends on (t) which gives the increment between time t-dt and t and which will be used as the stochastic process for the simulation (function)
%% OUTPUT:
%%          Physics2D: Structure containing variables concerning the physics of the problem in 2D (structure) (see Physics2D_Var2d.m)

function [Physics2D] = StochasticDispersion_Var2d(Method, Physics2D, StochasticDispersion, G, StochasticProcess)
%% Initializing the default potential
Default_Dispersion = @(W,FFTX,FFTY) -W.*FFTX.^2+FFTY.^2;

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
% IF there are 4 inputs and the potential is not defined
if (nargin >= 4) && (iscell(StochasticDispersion) == 0) && (isempty(StochasticDispersion) == 1) && (isempty(G) == 0)
    % FOR each component
    for n = 1:Method.Ncomponents
        Dispersion_function_index = []; % Initializing the temporary dispersion index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero
                if (G(n,m)~= 0)
                    % IF it is a diagonal term
                    if(n == m)
                    Physics2D.StochasticDispersion_function{n,m} = @(W,FFTX,FFTY)G(n,m)*Default_Dispersion(W,FFTX,FFTY); % Storing the dispersion as the correspondant value of G at the index added by the default dispersion
                    Dispersion_function_index = [Dispersion_function_index,m]; % Add the 'm' index in the temporary dispersion index
                    elseif (n ~= m)
                    Physics2D.StochasticDispersion_function{n,m} = @(W,FFTX,FFTY) 0; % Storing zero in the dispersion function matrix
                    end
                % ELSEIF the value of the matrix at the index is zero
                elseif (G(n,m) == 0)
                    Physics2D.StochasticDispersion_function{n,m} = @(W,FFTX,FFTY) 0; % Storing zero in the dispersion function matrix
                end
        end
        Physics2D.StochasticDispersion_function_Index{n} = Dispersion_function_index; % Store the dispersion index for the 'm' index
    end
% ELSEIF there are 4 inputs and the dispersion is defined but not a
% cell array
elseif (nargin >= 4) && (iscell(StochasticDispersion) == 0) && (isempty(StochasticDispersion) == 0) && (isempty(G) == 0)
    % FOR each component
    for n = 1:Method.Ncomponents
        Dispersion_function_index = []; % Initializing the temporary dispersion index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero
                if (G(n,m)~= 0)
                    Physics2D.StochasticDispersion_function{n,m} = @(W,FFTX,FFTY)G(n,m)*StochasticDispersion(W,FFTX,FFTY); % Storing the dispersion as the correspondant value of G at the index multiplied by the defined dispersion
                    Dispersion_function_index = [Dispersion_function_index,m]; % Add the 'm' index in the temporary dispersion index
                % ELSEIF the value of the matrix at the index is zero
                elseif (G(n,m) == 0)
                    Physics2D.StochasticDisersion_function{n,m} = @(W,FFTX,FFTY)0; % Storing zero in the dispersion function matrix
                end
        end
        Physics2D.StochasticDispersion_function_Index{n} = Dispersion_function_index; % Store the dispersion index for the 'm' index
    end
% ELSEIF there are 4 inputs and the dispersion is a cell array
elseif (nargin >= 4) && (iscell(StochasticDispersion) == 1) && (isempty(G) == 0)
    % FOR each component
    for n = 1:Method.Ncomponents
        Dispersion_function_index = []; % Initializing the temporary dispersion index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero and
                % the dispersion function is defined
                if (G(n,m)~= 0) && (isempty(StochasticDispersion{n,m}) == 0)
                    Physics2D.StochasticDispersion_function{n,m} = @(W,FFTX,FFTY)G(n,m)*StochasticDispersion{n,m}(W,FFTX,FFTY); % Storing the dispersion as the correspondant value of G at the index multiplied by the defined dispersion
                    Dispersion_function_index = [Dispersion_function_index,m]; % Add the 'm' index in the temporary dispersion index
                % ELSEIF the value of the matrix at the index is zero or
                % the dispersion function is not defined
                else
                    Physics2D.StochasticDispersion_function{n,m} = @(W,FFTX,FFTY)0; % Storing zero in the dispersion function matrix
                end
        end
        Physics2D.StochasticDispersion_function_Index{n} = Dispersion_function_index; % Store the dispersion index for the 'm' index
    end
    
% ELSEIF there are 4 inputs and the dispersion is not defined
elseif (nargin >= 4) && (iscell(StochasticDispersion) == 0) && (isempty(StochasticDispersion) == 1) && (isempty(G) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Dispersion_function_index = []; % Initializing the temporary dispersion index
        % FOR each component
        for m = 1:Method.Ncomponents
                    % IF it is a diagonal term
                    if(n == m)
                    Physics2D.StochasticDispersion_function{n,m} = @(W,FFTX,FFTY)Default_Dispersion(W,FFTX,FFTY); % Storing the dispersion as the correspondant value of G at the index added by the default dispersion
                    Dispersion_function_index = [Dispersion_function_index,m]; % Add the 'm' index in the temporary dispersion index
                    % ELSE if it is an extradiagonal term
                    else
                    Physics2D.StochasticDispersion_function{n,m} = @(W,FFTX,FFTY) 0; % Storing zero in the dispersion function matrix
                    end
        end
        Physics2D.StochasticDispersion_function_Index{n} = Dispersion_function_index; % Store the dispersion index for the 'm' index
    end
% ELSEIF there are 4 inputs and the dispersion is defined but not a
% cell array
elseif (nargin >= 4) && (iscell(StochasticDispersion) == 0) && (isempty(StochasticDispersion) == 0) && (isempty(G) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Dispersion_function_index = []; % Initializing the temporary dispersion index
        % FOR each component
        for m = 1:Method.Ncomponents
                    % IF it is a diagonal term
                    if (n == m)
                    Physics2D.StochasticDispersion_function{n,m} = @(W,FFTX,FFTY)StochasticDispersion(W,FFTX,FFTY); % Storing the dispersion as the correspondant value of G at the index multiplied by the defined dispersion
                    Dispersion_function_index = [Dispersion_function_index,m]; % Add the 'm' index in the temporary dispersion index
                    % ELSE if it is an extradiagonal term
                    else
                    Physics2D.StochasticDispersion_function{n,m} = @(W,FFTX,FFTY)0; % Storing zero in the dispersion function matrix
                    end
        end
        Physics2D.StochasticDispersion_function_Index{n} = Dispersion_function_index; % Store the dispersion index for the 'm' index
    end
% ELSEIF there are 4 inputs and the dispersion is a cell array
elseif (nargin >= 4) && (iscell(StochasticDispersion) == 1) && (isempty(G) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Dispersion_function_index = []; % Initializing the temporary dispersion index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the value of the matrix at the index is not zero and
                % the dispersion function is defined
                if (isempty(StochasticDispersion{n,m}) == 0)
                    Physics2D.StochasticDispersion_function{n,m} = @(W,FFTX,FFTY)StochasticDispersion{n,m}(W,FFTX,FFTY); % Storing the dispersion as the correspondant value of G at the index multiplied by the defined dispersion
                    Dispersion_function_index = [Dispersion_function_index,m]; % Add the 'm' index in the temporary dispersion index
                % ELSEIF the value of the matrix at the index is zero or
                % the dispersion function is not defined
                else
                    Physics2D.StochasticDispersion_function{n,m} = @(W,FFTX,FFTY)0; % Storing zero in the dispersion function matrix
                end
        end
        Physics2D.StochasticDispersion_function_Index{n} = Dispersion_function_index; % Store the dispersion index for the 'm' index
    end
% ELSEIF there are 3 inputs and the dispersion is not defined
elseif (nargin == 3) && (iscell(StochasticDispersion) == 0) && (isempty(StochasticDispersion) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Dispersion_function_index = []; % Initializing the temporary dispersion index
        % FOR each component
        for m = 1:Method.Ncomponents
            % IF it is a diagonal term
            if(n == m)
                Physics2D.StochasticDispersion_function{n,m} = @(W,FFTX,FFTY)G(n,m)*Default_Dispersion(W,FFTX,FFTY); % Storing the dispersion as the correspondant value of G at the index added by the default dispersion
                Dispersion_function_index = [Dispersion_function_index,m]; % Add the 'm' index in the temporary dispersion index
            elseif (n ~= m)
                Physics2D.StochasticDispersion_function{n,m} = @(W,FFTX,FFTY) 0; % Storing zero in the dispersion function matrix
            end
        end
        Physics2D.StochasticDispersion_function_Index{n} = Dispersion_function_index; % Store the dispersion index for the 'm' index
    end
% ELSEIF there are 3 inputs and the dispersion is defined but not a
% cell array
elseif (nargin == 3) && (iscell(StochasticDispersion) == 0) && (isempty(StochasticDispersion) == 0)
    % FOR each component
    for n = 1:Method.Ncomponents
        Dispersion_function_index = []; % Initializing the temporary dispersion index
        % FOR each component
        for m = 1:Method.Ncomponents
            Physics2D.StochasticDispersion_function{n,m} = @(W,FFTX,FFTY)StochasticDispersion(W,FFTX,FFTY); % Storing the dispersion as the correspondant value of G at the index multiplied by the defined dispersion
            Dispersion_function_index = [Dispersion_function_index,m]; % Add the 'm' index in the temporary dispersion index
        end
        Physics2D.StochasticDispersion_function_Index{n} = Dispersion_function_index; % Store the dispersion index for the 'm' index
    end
% ELSEIF there are 3 inputs and the dispersion is a cell array
elseif (nargin == 3) && (iscell(StochasticDispersion) == 1)
    % FOR each component
    for n = 1:Method.Ncomponents
        Dispersion_function_index = []; % Initializing the temporary dispersion index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF the dispersion function is defined
                if (isempty(StochasticDispersion{n,m}) == 0)
                    Physics2D.StochasticDispersion_function{n,m} = @(W,FFTX,FFTY)StochasticDispersion{n,m}(W,FFTX,FFTY); % Storing the dispersion as the correspondant value of G at the index multiplied by the defined dispersion
                    Dispersion_function_index = [Dispersion_function_index,m]; % Add the 'm' index in the temporary dispersion index
                % ELSEIF the dispersion function is not defined
                else
                    Physics2D.StochasticDispersion_function{n,m} = @(W,FFTX,FFTY)0; % Storing zero in the dispersion function matrix
                end
        end
        Physics2D.StochasticDispersion_function_Index{n} = Dispersion_function_index; % Store the dispersion index for the 'm' index
    end
% ELSEIF there are 2 inputs
elseif (nargin == 2)
    % FOR each component
    for n = 1:Method.Ncomponents
        Dispersion_function_index = []; % Initializing the temporary dispersion index
        % FOR each component
        for m = 1:Method.Ncomponents
                % IF it is a diagonal term
                if (n == m)
                    Physics2D.StochasticDispersion_function{n,m} = @(W,FFTX,FFTY)Default_Dispersion(W,FFTX,FFTY); % Storing the dispersion as the default potential
                    Dispersion_function_index = [Dispersion_function_index,m]; % Add the 'm' index in the temporary dispersion index
                % ELSEIF it is an extradiagonal term
                elseif (n ~= m)
                    Physics2D.StochasticDispersion_function{n,m} = @(W,FFTX,FFTY)0; % Storing zero in the dispersion function matrix
                end
        end
        Physics2D.StochasticDispersion_function_Index{n} = Dispersion_function_index; % Store the dispersion index for the 'm' index
    end
end