%% Creation of the 2D physics structure for the FFT
%% INPUTS:
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          Physics2D: Structure containing variables concerning the physics of the problem in 2D (structure) (see Physics2D_Var2d.m)
%%          FFTGeometry2D: Structure containing variables concerning the 2D geometry for the FFT (structure) (see FFTGeometry2D_Var2d.m)
%% OUTPUT:
%%          FFTPhysics2D: Structure containing variables concerning the physics of the problem in 2D in the FFT context (structure)

function [FFTPhysics2D] = FFTPhysics2D_Var2d(FFTPhi, Method, Physics2D, FFTGeometry2D, FFTOperators2D)
%% Initialization
FFTPhysics2D = Physics2D; % Copying the 2D physics

%% Computing the potential and gradients functions in the 2D geometry for the FFT
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        FFTPhysics2D.Potential{n,m} = Physics2D.Potential_function{n,m}(FFTGeometry2D.X,FFTGeometry2D.Y); % Computing the potential function in the 2D geometry for the FFT
        FFTPhysics2D.Dispersion{n,m} = Physics2D.Dispersion_function{n,m}(-1i*FFTOperators2D.Gx,-1i*FFTOperators2D.Gy); % Computing the dispersion function in the 2D geometry for the FFT
        FFTPhysics2D.TimeDispersion{n,m} = Physics2D.TimeDispersion_function{n,m}(Method.Iterations*Method.Deltat,-1i*FFTOperators2D.Gx,-1i*FFTOperators2D.Gy); % Computing the time-dependent dispersion function in the 2D geometry for the FFT
        if strcmp(Method.Type,'Splitting')
            FFTPhysics2D.Dispersionx{n,m} = Physics2D.Dispersion_function{n,m}(-1i*FFTOperators2D.Gx,0); % Computing the dispersion function in the 2D geometry for the FFT
            FFTPhysics2D.Dispersiony{n,m} = Physics2D.Dispersion_function{n,m}(0,-1i*FFTOperators2D.Gy); % Computing the dispersion function in the 2D geometry for the FFT
        end
        FFTPhysics2D.Nonlinearity{n,m} = FFTPhysics2D.Nonlinearity_function{n,m}(FFTPhi,FFTGeometry2D.X,FFTGeometry2D.Y); % Computing and storing the coupled nonlinearities between components
        FFTPhysics2D.FFTNonlinearity{n,m} = FFTPhysics2D.FFTNonlinearity_function{n,m}(FFTPhi,FFTGeometry2D.X,FFTGeometry2D.Y,-1i*FFTOperators2D.Gx,-1i*FFTOperators2D.Gy); % Computing and storing the coupled non-local nonlinearities between components
        FFTPhysics2D.GradientNLx{n,m} = FFTPhysics2D.GradientNLx_function{n,m}(FFTPhi,FFTGeometry2D.X,FFTGeometry2D.Y,-1i*FFTOperators2D.Gx,-1i*FFTOperators2D.Gy); % Computing and storing the coupled non-local gradient nonlinearities between components
        FFTPhysics2D.GradientNLy{n,m} = FFTPhysics2D.GradientNLy_function{n,m}(FFTPhi,FFTGeometry2D.X,FFTGeometry2D.Y,-1i*FFTOperators2D.Gx,-1i*FFTOperators2D.Gy); % Computing and storing the coupled non-local gradient nonlinearities between components
        FFTPhysics2D.TimeGradientx{n,m} = FFTPhysics2D.GradientNLx_function{n,m}(Method.Iterations*Method.Deltat,FFTGeometry2D.X,FFTGeometry2D.Y); % Computing and storing the coupled time-dependent gradient between components
        FFTPhysics2D.TimeGradienty{n,m} = FFTPhysics2D.GradientNLy_function{n,m}(Method.Iterations*Method.Deltat,FFTGeometry2D.X,FFTGeometry2D.Y); % Computing and storing the coupled time-dependent gradient between components
        FFTPhysics2D.TimePotential{n,m} = FFTPhysics2D.TimePotential_function{n,m}(Method.Iterations*Method.Deltat,FFTGeometry2D.X,FFTGeometry2D.Y); % Computing and storing the coupled time-dependent potential between components
        FFTPhysics2D.StochasticPotential{n,m} = Physics2D.StochasticPotential_function{n,m}(Physics2D.StochasticProcess_function(Method.Iterations*Method.Deltat,FFTGeometry2D.X,FFTGeometry2D.Y),FFTGeometry2D.X,FFTGeometry2D.Y) ; % Initializing the stochastic potential function
        FFTPhysics2D.StochasticDispersion{n,m} = Physics2D.StochasticDispersion_function{n,m}(Physics2D.StochasticProcess_function(Method.Iterations*Method.Deltat,FFTGeometry2D.X,FFTGeometry2D.Y),-1i*FFTOperators2D.Gx,-1i*FFTOperators2D.Gy) ; % Initializing the stochastic dispersion function
        FFTPhysics2D.TimePotentialExp{n,m} = FFTPhysics2D.TimePotential_function{n,m}(Method.Iterations*Method.Deltat,FFTGeometry2D.X,FFTGeometry2D.Y); % Computing and storing the coupled time-dependent potential between components
        FFTPhysics2D.TimePotentialImp{n,m} = FFTPhysics2D.TimePotential_function{n,m}((Method.Iterations + 1)*Method.Deltat,FFTGeometry2D.X,FFTGeometry2D.Y); % Computing and storing the coupled time-dependent potential between components
        FFTPhysics2D.TimeDispersionExp{n,m} = FFTPhysics2D.TimeDispersion_function{n,m}(Method.Iterations*Method.Deltat,-1i*FFTOperators2D.Gx,-1i*FFTOperators2D.Gy); % Computing and storing the coupled time-dependent dispersion between components
        FFTPhysics2D.TimeDispersionImp{n,m} = FFTPhysics2D.TimeDispersion_function{n,m}((Method.Iterations + 1)*Method.Deltat,-1i*FFTOperators2D.Gx,-1i*FFTOperators2D.Gy); % Computing and storing the coupled time-dependent dispersion between components
    end
end

% Computing the gradient in the x direction functions
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component where the gradient function is not null
    for m = Physics2D.Gradientx_function_Index{n}
        FFTPhysics2D.Gradientx{n,m} = Physics2D.Gradientx_function{n,m}(FFTGeometry2D.X,FFTGeometry2D.Y); % Computing the gradient in the x direction function in the 2D geometry for the FFT
    end
end

% Computing the gradient in the y direction functions
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component where the gradient function is not null
    for m = Physics2D.Gradienty_function_Index{n}
        FFTPhysics2D.Gradienty{n,m} = Physics2D.Gradienty_function{n,m}(FFTGeometry2D.X,FFTGeometry2D.Y); % Computing the gradient in the y direction function in the 2D geometry for the FFT
    end
end

% IF one wants to use the full Laplace preconditionner
if (strcmp(Method.Precond,'FLaplace') == 1)
    if (strcmp(Method.Type,'BESP') == 1)
        FFTPhysics2D.FPLaplace = BESPFPLaplace2d(Method, FFTPhysics2D, FFTGeometry2D); % Computing the full Laplace preconditioner for BESP
    elseif (strcmp(Method.Type,'CNSP') == 1)
        FFTPhysics2D.FPLaplace = CNSPFPLaplace2d(Method, FFTPhysics2D, FFTGeometry2D); % Computing the full Laplace preconditioner for CNSP
    elseif (strcmp(Method.Type,'Relaxation') == 1)
        FFTPhysics2D.FPLaplace = RSPFPLaplace2d(Method, FFTPhysics2D, FFTGeometry2D); % Computing the full Laplace preconditioner for BESP
    end
end

if strcmp(Method.Type,'Splitting')
    %% Initializing the exponantial operator
    MatcDiffx = cell(Method.Ncomponents); % Initializing the variable that will contain the differential operators in the x direction
    MatcDiffy = cell(Method.Ncomponents); % Initializing the variable that will contain the differential operators in the y direction

    % FOR each component
    for n = 1:Method.Ncomponents
        % FOR each component
        for m = 1:Method.Ncomponents
            MatcDiffx{n,m} = 0;
            MatcDiffy{n,m} = 0;
        end
    end

%% Constructing the cell matrix using the gradient operators
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = FFTPhysics2D.Gradientx_function_Index{n}
        MatcDiffx{n,m} = MatcDiffx{n,m} - 1i*Method.Deltat*FFTPhysics2D.Gradientx{n,m}.*FFTOperators2D.Gx; % Computing and adding the gradient operator
    end
    % FOR each component
    for m = FFTPhysics2D.Gradienty_function_Index{n}
        MatcDiffy{n,m} = MatcDiffy{n,m} - 1i*Method.Deltat*FFTPhysics2D.Gradienty{n,m}.*FFTOperators2D.Gy; % Computing and adding the gradient operator
    end
    % FOR each component where the dipersion is non null
    for m = FFTPhysics2D.Dispersion_function_Index{n}
        MatcDiffx{n,m} = MatcDiffx{n,m} - 1i*Method.Deltat*FFTPhysics2D.Dispersionx{n,m};
        MatcDiffy{n,m} = MatcDiffy{n,m} - 1i*Method.Deltat*FFTPhysics2D.Dispersiony{n,m};
    end
end

if (strcmp(Method.Splitting, 'Lie'))
    %% Lie splitting scheme on a single time step
    FFTPhysics2D.ExpMatcDiffx = expc(MatcDiffx); % Computing the exponential operator for the differential operators in the x direction
    FFTPhysics2D.ExpMatcDiffy = expc(MatcDiffy); % Computing the exponential operator for the differential operators in the y direction
elseif (strcmp(Method.Splitting, 'Strang'))
    %% Strang splitting scheme on a single time step
    FFTPhysics2D.ExpMatcDiffx = expc(Mulc(MatcDiffx,1/2)); % Computing the exponential operator for the differential operators in the x direction
    FFTPhysics2D.ExpMatcDiffy = expc(Mulc(MatcDiffy,1/2)); % Computing the exponential operator for the differential operators in the y direction
elseif (strcmp(Method.Splitting, 'Fourth'))
    %% Fourth-order splitting scheme on a single time step
    % Computing time-weights
    w_1 = 0.33780179798991440851;
    w_3 = -0.08780179798991440851;
    % Computing exponential operators
    FFTPhysics2D.ExpMatcDiffx1 = expc(Mulc(MatcDiffx,2*w_1)); % Computing the exponential operator for the differential operators in the x direction
    FFTPhysics2D.ExpMatcDiffx2 = expc(Mulc(MatcDiffx,2*w_3)); % Computing the exponential operator for the differential operators in the x direction
    FFTPhysics2D.ExpMatcDiffy1 = expc(Mulc(MatcDiffy,2*w_1)); % Computing the exponential operator for the differential operators in the y direction
    FFTPhysics2D.ExpMatcDiffy2 = expc(Mulc(MatcDiffy,2*w_3)); % Computing the exponential operator for the differential operators in the y direction
end
end