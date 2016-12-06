%% Computation of a step in time using the partial BESP-CNFG scheme
%% INPUTS:
%%          Phi: Initial wave functions in the 1D geometry for the FFT (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var1d.m)
%%          FFTGeometry1D: Structure containing variables concerning the geometry of the problem in 1D in the FFT context (structure) (see FFTGeometry1D_Var1d.m)
%%          FFTPhysics1D: Structure containing variables concerning the physics of the problem in 1D in the FFT context (structure) (see FFTPhysics1D_Var1d.m)
%%          FFTOperators1D: Structure containing the derivative FFT operators (structure) (see FFTOperators1D_Var1d.m)
%% OUTPUTS:
%%          Phi: Wave functions computated with the BESP method on a single step (cell array)
%%          iter: Total number of iterations before convergence (double)
%% FUNCTIONS USED:
%%           operator_Partial_BESP1d: To compute the wave functions by an iterative method (line 26)

function [Phi,iter] = Local_Partial_BESP_solution1d(Phi, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D)

%% Computing the explicit nonlinearity and the time-dependent potential
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component where the nonlinearity is non null
    for m = FFTPhysics1D.Nonlinearity_function_Index{n}
        FFTPhysics1D.Nonlinearity{n,m} = FFTPhysics1D.Nonlinearity_function{n,m}(Phi,FFTGeometry1D.X); % Computing and storing the coupled nonlinearities between components
    end
    % FOR each component where the non-local nonlinearity is non null
    for m = FFTPhysics1D.FFTNonlinearity_function_Index{n}
        FFTPhysics1D.FFTNonlinearity{n,m} = FFTPhysics1D.FFTNonlinearity_function{n,m}(Phi,FFTGeometry1D.X,-1i*FFTOperators1D.Gx); % Computing and storing the coupled non-local nonlinearities between components
    end
    % FOR each component where the time-dependent potential is non null
    for m = FFTPhysics1D.TimePotential_function_Index{n}
        FFTPhysics1D.TimePotential{n,m} = FFTPhysics1D.TimePotential_function{n,m}(Method.Iterations*Method.Deltat,FFTGeometry1D.X); % Computing and storing the coupled time-dependent potential between components
    end
end

%% Computing the wave functions after a single time step using the direct iterative method for BESP-CNFG
[Phi,iter] = operator_Partial_BESP1d(Phi, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D); %Using the direct iterative method to solve a single step in time of the BESP-CNFG scheme