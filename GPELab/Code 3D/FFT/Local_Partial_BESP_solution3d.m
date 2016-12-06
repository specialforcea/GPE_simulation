%% Computation of a step in time using the partial BESP-CNFG scheme
%% INPUTS:
%%          Phi: Initial wave functions in the 3D geometry for the FFT (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var3d.m)
%%          FFTGeometry3D: Structure containing variables concerning the geometry of the problem in 3D in the FFT context (structure) (see FFTGeometry3D_Var3d.m)
%%          FFTPhysics3D: Structure containing variables concerning the physics of the problem in 3D in the FFT context (structure) (see FFTPhysics3D_Var3d.m)
%%          FFTOperators3D: Structure containing the derivative FFT operators (structure) (see FFTOperators3D_Var3d.m)
%% OUTPUTS:
%%          Phi: Wave functions computated with the BESP method on a single step (cell array)
%%          iter: Total number of iterations before convergence (double)
%% FUNCTIONS USED:
%%           operator_Partial_BESP3d: To compute the wave functions by an iterative method (line 30)

function [Phi,iter] = Local_Partial_BESP_solution3d(Phi, Method, FFTGeometry3D, FFTPhysics3D, FFTOperators3D)

%% Computing the explicit nonlinearity
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component where the nonlinearity is non null
    for m = FFTPhysics3D.Nonlinearity_function_Index{n}
        FFTPhysics3D.Nonlinearity{n,m} = FFTPhysics3D.Nonlinearity_function{n,m}(Phi,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z); % Computing and storing the coupled nonlinearities between components
    end
    % FOR each component where the non-local nonlinearity is non null
    for m = FFTPhysics3D.FFTNonlinearity_function_Index{n}
        FFTPhysics3D.FFTNonlinearity{n,m} = FFTPhysics3D.FFTNonlinearity_function{n,m}(Phi,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z,-1i*FFTOperators3D.Gx,-1i*FFTOperators3D.Gy,-1i*FFTOperators3D.Gz); % Computing and storing the coupled non-local nonlinearities between components
    end
    % FOR each component where the time-dependent potential is non null
    for m = FFTPhysics3D.TimePotential_function_Index{n}
        FFTPhysics3D.TimePotential{n,m} = FFTPhysics3D.TimePotential_function{n,m}(Method.Iterations*Method.Deltat,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z); % Computing and storing the coupled time-dependent potentials between components
    end
end

%% Computing the wave functions after a single time step using the direct iterative method for BESP-CNFG
[Phi,iter] = operator_Partial_BESP3d(Phi, Method, FFTGeometry3D, FFTPhysics3D, FFTOperators3D); %Using the direct iterative method to solve a single step in time of the BESP-CNFG scheme