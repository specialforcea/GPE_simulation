%% Creation of the 3D physics structure for the FFT
%% INPUTS:
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var3d.m)
%%          Physics3D: Structure containing variables concerning the physics of the problem in 3D (structure) (see Physics3D_Var3d.m)
%%          FFTGeometry3D: Structure containing variables concerning the 3D geometry for the FFT (structure) (see FFTGeometry3D_Var3d.m)
%% OUTPUT:
%%          FFTPhysics3D: Structure containing variables concerning the physics of the problem in 3D in the FFT context (structure)

function [FFTPhysics3D] = FFTPhysics3D_Var3d(Phi, Method, Physics3D, FFTGeometry3D, FFTOperators3D)
%% Initialization
FFTPhysics3D = Physics3D; % Copying the 3D physics

%% Computing the potential and gradients functions in the 3D geometry for the FFT
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        FFTPhysics3D.Potential{n,m} = Physics3D.Potential_function{n,m}(FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z); % Computing the potential function in the 3D geometry for the FFT
        FFTPhysics3D.Dispersion{n,m} = Physics3D.Dispersion_function{n,m}(-1i*FFTOperators3D.Gx,-1i*FFTOperators3D.Gy,-1i*FFTOperators3D.Gz); % Computing the dispersion function in the 3D geometry for the FFT
        if strcmp(Method.Type,'Splitting')
        FFTPhysics3D.Dispersionx{n,m} = Physics3D.Dispersion_function{n,m}(-1i*FFTOperators3D.Gx,0,0); % Computing the dispersion function in the 3D geometry for the FFT
        FFTPhysics3D.Dispersiony{n,m} = Physics3D.Dispersion_function{n,m}(0,-1i*FFTOperators3D.Gy,0); % Computing the dispersion function in the 3D geometry for the FFT
        FFTPhysics3D.Dispersionz{n,m} = Physics3D.Dispersion_function{n,m}(0,0,-1i*FFTOperators3D.Gz); % Computing the dispersion function in the 3D geometry for the FFT      
        end
        FFTPhysics3D.TimeDispersion{n,m} = Physics3D.TimeDispersion_function{n,m}(Method.Deltat*Method.Iterations, FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z); % Computing the time dispersion function in the 3D geometry for the FFT
        FFTPhysics3D.TimeDispersionImp{n,m} = Physics3D.TimeDispersion_function{n,m}(Method.Deltat*Method.Iterations, FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z); % Computing the time dispersion function in the 3D geometry for the FFT
        FFTPhysics3D.TimeDispersionExp{n,m} = Physics3D.TimeDispersion_function{n,m}(Method.Deltat*Method.Iterations, FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z); % Computing the time dispersion function in the 3D geometry for the FFT
        FFTPhysics3D.TimePotential{n,m} = Physics3D.TimePotential_function{n,m}(Method.Deltat*Method.Iterations, FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z); % Computing the potential function in the 3D geometry for the FFT
        FFTPhysics3D.Nonlinearity{n,m} = Physics3D.Nonlinearity_function{n,m}(Phi, FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z); % Computing the potential function in the 3D geometry for the FFT
        FFTPhysics3D.FFTNonlinearity{n,m} = Physics3D.FFTNonlinearity_function{n,m}(Phi, FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z,-1i*FFTOperators3D.Gx,-1i*FFTOperators3D.Gy,-1i*FFTOperators3D.Gz); % Computing the potential function in the 3D geometry for the FFT
        FFTPhysics3D.StochasticPotential{n,m} = Physics3D.StochasticPotential_function{n,m}(Physics3D.StochasticProcess_function(Method.Iterations*Method.Deltat,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z),FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z) ; % Initializing the stochastic potential function
        FFTPhysics3D.StochasticDispersion{n,m} = Physics3D.StochasticDispersion_function{n,m}(Physics3D.StochasticProcess_function(Method.Iterations*Method.Deltat,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z),-1i*FFTOperators3D.Gx,-1i*FFTOperators3D.Gy,-1i*FFTOperators3D.Gz) ; % Initializing the stochastic dispersion function
        FFTPhysics3D.TimePotentialExp{n,m} = FFTPhysics3D.TimePotential_function{n,m}(Method.Iterations*Method.Deltat,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z); % Computing and storing the coupled time-dependent potential between components
        FFTPhysics3D.TimePotentialImp{n,m} = FFTPhysics3D.TimePotential_function{n,m}((Method.Iterations + 1)*Method.Deltat,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z); % Computing and storing the coupled time-dependent potential between components
    end
end

% Computing the gradient in the x direction functions
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component where the potential function is not null
    for m = Physics3D.Gradientx_function_Index{n}
        FFTPhysics3D.Gradientx{n,m} = Physics3D.Gradientx_function{n,m}(FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z); % Computing the gradient in the x direction function in the 3D geometry for the FFT
    end
end

% Computing the gradient in the y direction functions
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component where the potential function is not null
    for m = Physics3D.Gradienty_function_Index{n}
        FFTPhysics3D.Gradienty{n,m} = Physics3D.Gradienty_function{n,m}(FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z); % Computing the gradient in the y direction function in the 3D geometry for the FFT
    end
end

% Computing the gradient in the z direction functions
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component where the potential function is not null
    for m = Physics3D.Gradientz_function_Index{n}
        FFTPhysics3D.Gradientz{n,m} = Physics3D.Gradientz_function{n,m}(FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z); % Computing the gradient in the z direction function in the 3D geometry for the FFT
    end
end


% IF one wants to use the full Laplace preconditionner
if (strcmp(Method.Precond,'FLaplace') == 1)
    if (strcmp(Method.Type,'BESP') == 1)
        FFTPhysics3D.FPLaplace = BESPFPLaplace3d(Method, FFTPhysics3D, FFTGeometry3D); % Computing the full Laplace preconditioner for BESP
    elseif (strcmp(Method.Type,'CNSP') == 1)
        FFTPhysics3D.FPLaplace = CNSPFPLaplace3d(Method, FFTPhysics3D, FFTGeometry3D); % Computing the full Laplace preconditioner for CNSP
    elseif (strcmp(Method.Type,'Relaxation') == 1)
        FFTPhysics3D.FPLaplace = RSPFPLaplace3d(Method, FFTPhysics3D, FFTGeometry3D); % Computing the full Laplace preconditioner for RSP
    end
end


if strcmp(Method.Type,'Splitting')
MatcDiffx = cell(Method.Ncomponents); % Initializing the variable that will contain the differential operators in the x direction
MatcDiffy = cell(Method.Ncomponents); % Initializing the variable that will contain the differential operators in the y direction
MatcDiffz = cell(Method.Ncomponents); % Initializing the variable that will contain the differential operators in the z direction
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        MatcDiffx{n,m} = 0;
        MatcDiffy{n,m} = 0;
        MatcDiffz{n,m} = 0;
    end
end

%% Constructing the cell matrix using the potential and nonlinear operators
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = FFTPhysics3D.Gradientx_function_Index{n}
        MatcDiffx{n,m} = MatcDiffx{n,m} - 1i*Method.Deltat*FFTPhysics3D.Gradientx{n,m}.*FFTOperators3D.Gx; % Computing and adding the potential operator
    end
    % FOR each component
    for m = FFTPhysics3D.Gradienty_function_Index{n}
        MatcDiffy{n,m} = MatcDiffy{n,m} - 1i*Method.Deltat*FFTPhysics3D.Gradienty{n,m}.*FFTOperators3D.Gy; % Computing and adding the potential operator
    end
    % FOR each component
    for m = FFTPhysics3D.Gradientz_function_Index{n}
        MatcDiffz{n,m} = MatcDiffz{n,m} - 1i*Method.Deltat*FFTPhysics3D.Gradientz{n,m}.*FFTOperators3D.Gz; % Computing and adding the potential operator
    end
    % FOR each component where the dipersion is non null
    for m = FFTPhysics3D.Dispersion_function_Index{n}
        MatcDiffx{n,m} = MatcDiffx{n,m} - 1i*Method.Deltat*FFTPhysics3D.Dispersionx{n,m};
        MatcDiffy{n,m} = MatcDiffy{n,m} - 1i*Method.Deltat*FFTPhysics3D.Dispersiony{n,m};
        MatcDiffz{n,m} = MatcDiffz{n,m} - 1i*Method.Deltat*FFTPhysics3D.Dispersionz{n,m};
    end 
end

if (strcmp(Method.Splitting, 'Lie'))
%% Lie splitting scheme on a single time step
FFTPhysics3D.ExpMatcDiffx = expc(MatcDiffx); % Computing the exponential operator for the differential operators in the x direction
FFTPhysics3D.ExpMatcDiffy = expc(MatcDiffy); % Computing the exponential operator for the differential operators in the y direction
FFTPhysics3D.ExpMatcDiffz = expc(MatcDiffz); % Computing the exponential operator for the differential operators in the z direction

elseif (strcmp(Method.Splitting, 'Strang'))
%% Strang splitting scheme on a single time step
FFTPhysics3D.ExpMatcDiffx = expc(Mulc(MatcDiffx,1/2)); % Computing the exponential operator for the differential operators in the x direction
FFTPhysics3D.ExpMatcDiffy = expc(Mulc(MatcDiffy,1/2)); % Computing the exponential operator for the differential operators in the y direction
FFTPhysics3D.ExpMatcDiffz = expc(Mulc(MatcDiffz,1/2)); % Computing the exponential operator for the differential operators in the z direction

elseif (strcmp(Method.Splitting, 'Fourth'))
%% Fourth-order splitting scheme on a single time step
% Computing time-weights
w_1 = 0.33780179798991440851;
w_3 = -0.08780179798991440851;

% Computing exponential operators
FFTPhysics3D.ExpMatcDiffx1 = expc(Mulc(MatcDiffx,2*w_1)); % Computing the exponential operator for the differential operators in the x direction
FFTPhysics3D.ExpMatcDiffx2 = expc(Mulc(MatcDiffx,2*w_3)); % Computing the exponential operator for the differential operators in the x direction
FFTPhysics3D.ExpMatcDiffy1 = expc(Mulc(MatcDiffy,2*w_1)); % Computing the exponential operator for the differential operators in the y direction
FFTPhysics3D.ExpMatcDiffy2 = expc(Mulc(MatcDiffy,2*w_3)); % Computing the exponential operator for the differential operators in the y direction
FFTPhysics3D.ExpMatcDiffz1 = expc(Mulc(MatcDiffz,2*w_1)); % Computing the exponential operator for the differential operators in the z direction
FFTPhysics3D.ExpMatcDiffz2 = expc(Mulc(MatcDiffz,2*w_3)); % Computing the exponential operator for the differential operators in the z direction

end
end