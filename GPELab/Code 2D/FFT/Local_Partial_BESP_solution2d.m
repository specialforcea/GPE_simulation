%% Computation of a step in time using the partial BESP-CNFG scheme
%% INPUTS:
%%          Phi: Initial wave functions in the 2D geometry for the FFT (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          FFTGeometry2D: Structure containing variables concerning the geometry of the problem in 2D in the FFT context (structure) (see FFTGeometry2D_Var2d.m)
%%          FFTPhysics2D: Structure containing variables concerning the physics of the problem in 2D in the FFT context (structure) (see FFTPhysics2D_Var2d.m)
%%          FFTOperators2D: Structure containing the derivative FFT operators (structure) (see FFTOperators2D_Var2d.m)
%% OUTPUTS:
%%          Phi: Wave functions computated with the BESP method on a single step (cell array)
%%          iter: Total number of iterations before convergence (double)
%% FUNCTIONS USED:
%%           operator_Partial_BESP2d: To compute the wave functions by an iterative method (line 26)

function [Phi,iter] = Local_Partial_BESP_solution2d(Phi, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D)

%% Computing the explicit nonlinearity, non-local nonlinearity and the time-dependent potential
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component where the nonlinearity is non null
    for m = FFTPhysics2D.Nonlinearity_function_Index{n}
        FFTPhysics2D.Nonlinearity{n,m} = FFTPhysics2D.Nonlinearity_function{n,m}(Phi,FFTGeometry2D.X,FFTGeometry2D.Y); % Computing and storing the coupled nonlinearities between components
    end
    % FOR each component where the non-local nonlinearity is non null
    for m = FFTPhysics2D.FFTNonlinearity_function_Index{n}
        FFTPhysics2D.FFTNonlinearity{n,m} = FFTPhysics2D.FFTNonlinearity_function{n,m}(Phi,FFTGeometry2D.X,FFTGeometry2D.Y,-1i*FFTOperators2D.Gx,-1i*FFTOperators2D.Gy); % Computing and storing the coupled non-local nonlinearities between components
    end
    % FOR each component where the time-dependent potential is non null
    for m = FFTPhysics2D.TimePotential_function_Index{n}
        FFTPhysics2D.TimePotential{n,m} = FFTPhysics2D.TimePotential_function{n,m}(Method.Iterations*Method.Deltat,FFTGeometry2D.X,FFTGeometry2D.Y); % Computing and storing the coupled time-dependent potential between components
    end
end

%% Computing the wave functions after a single time step using the direct iterative method for BESP-CNFG
[Phi,iter] = operator_Partial_BESP2d(Phi, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D); %Using the direct iterative method to solve a single step in time of the BESP-CNFG scheme