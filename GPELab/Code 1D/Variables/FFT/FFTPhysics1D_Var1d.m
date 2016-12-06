%% Creation of the 1D physics structure for the FFT
%% INPUTS:
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var1d.m)
%%          Physics1D: Structure containing variables concerning the physics of the problem in 1D (structure) (see Physics1D_Var1d.m)
%%          FFTGeometry1D: Structure containing variables concerning the 1D geometry for the FFT (structure) (see FFTGeometry1D_Var1d.m)
%% OUTPUT:
%%          FFTPhysics1D: Structure containing variables concerning the physics of the problem in 1D in the FFT context (structure)

function [FFTPhysics1D] = FFTPhysics1D_Var1d(FFTPhi, Method, Physics1D, FFTGeometry1D, FFTOperators1D)
%% Initialization
FFTPhysics1D = Physics1D; % Copying the 1D physics

%% Computing the potential and gradients functions in the 1D geometry for the FFT, and setting the diagonal nonlinearity and time-dependent potential to zero in the case of a Thomas-Fermi preconditionner
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component where the potential function is not null
    for m = 1:Method.Ncomponents
        FFTPhysics1D.Potential{n,m} = Physics1D.Potential_function{n,m}(FFTGeometry1D.X); % Computing the potential function in the 1D geometry for the FFT
        FFTPhysics1D.Dispersion{n,m} = Physics1D.Dispersion_function{n,m}(-1i*FFTOperators1D.Gx); % Computing the dispersion function in the 1D geometry for the FFT
        FFTPhysics1D.TimeDispersion{n,m} = Physics1D.TimeDispersion_function{n,m}(Method.Iterations*Method.Deltat,-1i*FFTOperators1D.Gx); % Computing the time-dependent dispersion function in the 1D geometry for the FFT
        FFTPhysics1D.TimePotential{n,m} = Physics1D.TimePotential_function{n,m}(Method.Iterations*Method.Deltat,FFTGeometry1D.X); % Computing the time dependent potential function in the 1D geometry for the FFT
        FFTPhysics1D.StochasticPotential{n,m} = Physics1D.StochasticPotential_function{n,m}(Physics1D.StochasticProcess_function(Method.Iterations*Method.Deltat,FFTGeometry1D.X),FFTGeometry1D.X) ; % Initializing the stochastic potential function
        FFTPhysics1D.StochasticDispersion{n,m} = Physics1D.StochasticDispersion_function{n,m}(Physics1D.StochasticProcess_function(Method.Iterations*Method.Deltat,FFTGeometry1D.X),-1i*FFTOperators1D.Gx) ; % Initializing the stochastic potential function
        FFTPhysics1D.TimePotential{n,m} = Physics1D.TimePotential_function{n,m}(Method.Iterations*Method.Deltat,FFTGeometry1D.X); % Computing the time dependent potential function in the 1D geometry for the FFT
        FFTPhysics1D.Nonlinearity{n,m} = Physics1D.Nonlinearity_function{n,m}(FFTPhi,FFTGeometry1D.X); % Computing the nonlinearity in the 1D geometry for the FFT
        FFTPhysics1D.FFTNonlinearity{n,m} = Physics1D.FFTNonlinearity_function{n,m}(FFTPhi,FFTGeometry1D.X, -1i*FFTOperators1D.Gx); % Computing the nonlinearity in the 1D geometry for the FFT
        FFTPhysics1D.TimePotentialExp{n,m} = FFTPhysics1D.TimePotential_function{n,m}(Method.Iterations*Method.Deltat,FFTGeometry1D.X); % Computing and storing the coupled time-dependent potential between components
        FFTPhysics1D.TimePotentialImp{n,m} = FFTPhysics1D.TimePotential_function{n,m}((Method.Iterations + 1)*Method.Deltat,FFTGeometry1D.X); % Computing and storing the coupled time-dependent potential between components
    end
end

% Computing the gradient in the x direction functions
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component where the potential function is not null
    for m = Physics1D.Gradientx_function_Index{n}
        FFTPhysics1D.Gradientx{n,m} = Physics1D.Gradientx_function{n,m}(FFTGeometry1D.X); % Computing the gradient in the x direction function in the 1D geometry for the FFT
    end
end

% IF one wants to use the full Laplace preconditionner
if (strcmp(Method.Precond,'FLaplace') == 1)
    if (strcmp(Method.Type,'BESP') == 1)
        FFTPhysics1D.FPLaplace = BESPFPLaplace1d(Method, FFTPhysics1D, FFTGeometry1D); % Computing the full Laplace preconditioner for BESP
    elseif (strcmp(Method.Type,'CNSP') == 1)
        FFTPhysics1D.FPLaplace = CNSPFPLaplace1d(Method, FFTPhysics1D, FFTGeometry1D); % Computing the full Laplace preconditioner for BESP
    end
end

if strcmp(Method.Type,'Splitting')
MatcDiffx = cell(Method.Ncomponents); % Initializing the variable that will contain the differential operators
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        MatcDiffx{n,m} = 0;
    end
end

%% Constructing the cell matrix using the potential and nonlinear operators
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = FFTPhysics1D.Gradientx_function_Index{n}
        MatcDiffx{n,m} = MatcDiffx{n,m} - 1i*Method.Deltat*FFTPhysics1D.Gradientx{n,m}.*FFTOperators1D.Gx; % Computing and adding the gradient operator
    end
    % FOR each component where the stochastic dispersion is non null
    for m = FFTPhysics1D.Dispersion_function_Index{n}
        MatcDiffx{n,m} = MatcDiffx{n,m} - 1i*Method.Deltat*FFTPhysics1D.Dispersion{n,m}; % Computing and adding the dispersion operator
    end
end

if (strcmp(Method.Splitting, 'Lie'))
%% Lie splitting scheme on a single time step
FFTPhysics1D.ExpMatcDiffx = expc(MatcDiffx); % Computing the exponential operator for the differential operators in the x direction

elseif (strcmp(Method.Splitting, 'Strang'))
%% Strang splitting scheme on a single time step
FFTPhysics1D.ExpMatcDiffx = expc(Mulc(MatcDiffx,1/2)); % Computing the exponential operator for the differential operators in the x direction

elseif (strcmp(Method.Splitting, 'Fourth'))
%% Fourth-order splitting scheme on a single time step
% Computing time-weights
w_2 = 0.67560359597982881702;
w_4 = -0.85120719795965763405;

% Computing exponential operators
FFTPhysics1D.ExpMatcDiffx1 = expc(Mulc(MatcDiffx,2*w_2)); % Computing the exponential operator for the differential operators in the x direction
FFTPhysics1D.ExpMatcDiffx2 = expc(Mulc(MatcDiffx,2*w_4)); % Computing the exponential operator for the differential operators in the x direction

end
end