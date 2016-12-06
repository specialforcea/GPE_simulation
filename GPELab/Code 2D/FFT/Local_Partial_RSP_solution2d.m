%% Computation of a step in time using the partial RSP scheme
%% INPUT:
%%          Phi: Initial wave functions in the 2D geometry for the FFT (cell array)
%%          Psi: Relaxation variable corresponding to the nonlinearity (cell array)
%%          FFTPsi: Relaxation variable corresponding to the non local nonlinearity (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          FFTGeometry2D: Structure containing variables concerning the geometry of the problem in 2D in the FFT context (structure) (see FFTGeometry2D_Var2d.m)
%%          FFTPhysics2D: Structure containing variables concerning the physics of the problem in 2D in the FFT context (structure) (see FFTPhysics2D_Var2d.m)
%%          FFTOperators2D: Structure containing the derivative FFT operators (structure) (see FFTOperators2D_Var2d.m)
%% OUTPUTS:
%%          Phi: Wave functions computated with the BESP method on a single step (cell array)
%%          Psi: Relaxation variable corresponding to the nonlinearity (cell array)
%%          FFTPsi: Relaxation variable corresponding to the non local nonlinearity (cell array)
%%          [flag,relres,iter,resvec]: Outputs from bicgstab (vector) (see bicgstab.m)
%% FUNCTIONS USED:
%%          bicgstab: To compute the wave functions by an iterative method (line 53, 60 and 66)
%%          operator_Full_CNSP_PLaplacian2d: To compute the implicit part of the RSP scheme with the Laplace preconditioner already applied (line 53)
%%          rhs_Full_RSP_PLaplacian2d: To compute the explicit part of the RSP scheme with the Laplace preconditioner already applied (line 54)
%%          operator_Full_CNSP_PThomasFermi2d: To compute the implicit part of the RSP scheme with the Thomas-Fermi preconditioner already applied (line 60)
%%          rhs_Full_RSP_PThomasFermi2d: To compute the explicit part of the RSP scheme with the Thomas-Fermi preconditioner already applied (line 61)
%%          operator_Full_CNSP2d: To compute the implicit part of the RSP scheme (line 66)
%%          rhs_Full_RSP2d: To compute the explicit part of the RSP scheme (line 67)

function [Phi, Psi ,FFTPsi, iter] = Local_Partial_RSP_solution2d(Phi, Psi, FFTPsi, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D)

%% Computing the explicit nonlinearity, non-local nonlinearity and the time-depend potential
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component where the nonlinearity is non null
    for m = FFTPhysics2D.Nonlinearity_function_Index{n}
        FFTPhysics2D.Nonlinearity{n,m} = FFTPhysics2D.Nonlinearity_function{n,m}(Phi,FFTGeometry2D.X,FFTGeometry2D.Y); % Computing and storing the coupled nonlinearities between components
        Psi{n,m} = 2*FFTPhysics2D.Nonlinearity{n,m} - Psi{n,m}; % Computing the relaxation variable corresponding to the local nonlinearities
    end
    % FOR each component where the non-local nonlinearity is non null
    for m = FFTPhysics2D.FFTNonlinearity_function_Index{n}
        FFTPhysics2D.FFTNonlinearity{n,m} = FFTPhysics2D.FFTNonlinearity_function{n,m}(Phi,FFTGeometry2D.X,FFTGeometry2D.Y,-1i*FFTOperators2D.Gx,-1i*FFTOperators2D.Gy); % Computing and storing the coupled non-local nonlinearities between components
        FFTPsi{n,m} = 2*FFTPhysics2D.FFTNonlinearity{n,m} - FFTPsi{n,m}; % Computing the relaxation variable corresponding to the non-local nonlinearities
    end
    % FOR each component where the time-dependent potential is non null
    for m = FFTPhysics2D.TimePotential_function_Index{n}
        FFTPhysics2D.TimePotentialExp{n,m} = FFTPhysics2D.TimePotential_function{n,m}(Method.Iterations*Method.Deltat,FFTGeometry2D.X,FFTGeometry2D.Y); % Computing and storing the coupled time-dependent potential between components
        FFTPhysics2D.TimePotentialImp{n,m} = FFTPhysics2D.TimePotential_function{n,m}((Method.Iterations + 1)*Method.Deltat,FFTGeometry2D.X,FFTGeometry2D.Y); % Computing and storing the coupled time-dependent potential between components
    end
    % FOR each component where the stochastic dispersion is non null
    for m = FFTPhysics2D.StochasticDispersion_function_Index{n}
        if (iscell(FFTPhysics2D.StochasticProcess_function))
            for m_noise = 1:length(FFTPhysics2D.StochasticProcess_function)
                FFTPhysics2D.StochasticProcess{m_noise} = (FFTPhysics2D.StochasticProcess_function{m_noise}(Method.Iterations*Method.Deltat,FFTGeometry2D.X,FFTGeometry2D.Y)-FFTPhysics2D.StochasticProcess_function{m_noise}((Method.Iterations - 1)*Method.Deltat,FFTGeometry2D.X,FFTGeometry2D.Y))/Method.Deltat;
            end
        else
            FFTPhysics2D.StochasticProcess = (FFTPhysics2D.StochasticProcess_function(Method.Iterations*Method.Deltat,FFTGeometry2D.X,FFTGeometry2D.Y)-FFTPhysics2D.StochasticProcess_function((Method.Iterations - 1)*Method.Deltat,FFTGeometry2D.X,FFTGeometry2D.Y))/Method.Deltat;
        end
        FFTPhysics2D.StochasticDispersion{n,m} = FFTPhysics2D.StochasticDispersion_function{n,m}(FFTPhysics2D.StochasticProcess,-1i*FFTOperators2D.Gx,-1i*FFTOperators2D.Gy); % Computing and storing the stochastic potential
    end
end

%% Computing the wave functions after a single time step using Krylov iterative method for RSP
%If the computation is dynamic
if (strcmp(Method.Computation,'Dynamic'))
    Method.Deltat = 1i*Method.Deltat;
end

[Phi, Psi, FFTPsi, iter] = operator_Partial_RSP2d(Phi, Psi, FFTPsi, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D);