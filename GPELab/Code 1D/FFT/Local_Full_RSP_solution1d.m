%% Computation of a step in time using the RSP scheme
%% INPUT:
%%          Phi: Initial wave functions in the 2D geometry for the FFT (cell array)
%%          Psi: Relaxation variable corresponding to the nonlinearity (cell array)
%%          FFTPsi: Relaxation variable corresponding to the non local nonlinearity (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var1d.m)
%%          FFTGeometry1D: Structure containing variables concerning the geometry of the problem in 1D in the FFT context (structure) (see FFTGeometry1D_Var1d.m)
%%          FFTPhysics1D: Structure containing variables concerning the physics of the problem in 1D in the FFT context (structure) (see FFTPhysics1D_Var1d.m)
%%          FFTOperators1D: Structure containing the derivative FFT operators (structure) (see FFTOperators1D_Var1d.m)
%% OUTPUTS:
%%          Phi: Wave functions computated with the BESP method on a single step (cell array)
%%          Psi: Relaxation variable corresponding to the nonlinearity (cell array)
%%          FFTPsi: Relaxation variable corresponding to the non local nonlinearity (cell array)
%%          [flag,relres,iter,resvec]: Outputs from bicgstab (vector) (see bicgstab.m)
%% FUNCTIONS USED:
%%          bicgstab: To compute the wave functions by an iterative method (line 53, 60 and 66)
%%          operator_Full_CNSP_PLaplacian1d: To compute the implicit part of the RSP scheme with the Laplace preconditioner already applied (line 53)
%%          rhs_Full_RSP_PLaplacian1d: To compute the explicit part of the RSP scheme with the Laplace preconditioner already applied (line 54)
%%          operator_Full_CNSP_PThomasFermi1d: To compute the implicit part of the RSP scheme with the Thomas-Fermi preconditioner already applied (line 60)
%%          rhs_Full_RSP_PThomasFermi1d: To compute the explicit part of the RSP scheme with the Thomas-Fermi preconditioner already applied (line 61)
%%          operator_Full_CNSP1d: To compute the implicit part of the RSP scheme (line 66)
%%          rhs_Full_RSP1d: To compute the explicit part of the RSP scheme (line 67)

function [Phi, Psi ,FFTPsi, flag, relres, iter, resvec] = Local_Full_RSP_solution1d(Phi, Psi, FFTPsi, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D)

%% Computing the explicit nonlinearity, non-local nonlinearity and the time-depend potential
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component where the nonlinearity is non null
    for m = FFTPhysics1D.Nonlinearity_function_Index{n}
        FFTPhysics1D.Nonlinearity{n,m} = FFTPhysics1D.Nonlinearity_function{n,m}(Phi,FFTGeometry1D.X); % Computing and storing the coupled nonlinearities between components
        Psi{n,m} = 2*FFTPhysics1D.Nonlinearity{n,m} - Psi{n,m}; % Computing the relaxation variable corresponding to the local nonlinearities
    end
    % FOR each component where the non-local nonlinearity is non null
    for m = FFTPhysics1D.FFTNonlinearity_function_Index{n}
        FFTPhysics1D.FFTNonlinearity{n,m} = FFTPhysics1D.FFTNonlinearity_function{n,m}(Phi,FFTGeometry1D.X,-1i*FFTOperators1D.Gx); % Computing and storing the coupled non-local nonlinearities between components
        FFTPsi{n,m} = 2*FFTPhysics1D.FFTNonlinearity{n,m} - FFTPsi{n,m}; % Computing the relaxation variable corresponding to the non-local nonlinearities
    end
    % FOR each component where the time-dependent potential is non null
    for m = FFTPhysics1D.TimePotential_function_Index{n}
        FFTPhysics1D.TimePotentialExp{n,m} = FFTPhysics1D.TimePotential_function{n,m}((Method.Iterations-1)*Method.Deltat,FFTGeometry1D.X); % Computing and storing the coupled time-dependent potential between components
        FFTPhysics1D.TimePotentialImp{n,m} = FFTPhysics1D.TimePotential_function{n,m}(Method.Iterations*Method.Deltat,FFTGeometry1D.X); % Computing and storing the coupled time-dependent potential between components
    end
    % FOR each component where the time-dependent dispersion is non null
    for m = FFTPhysics1D.TimeDispersion_function_Index{n}
        FFTPhysics1D.TimeDispersionExp{n,m} = FFTPhysics1D.TimeDispersion_function{n,m}((Method.Iterations-1)*Method.Deltat,-1i*FFTOperators1D.Gx); % Computing and storing the coupled time-dependent potential between components
        FFTPhysics1D.TimeDispersionImp{n,m} = FFTPhysics1D.TimeDispersion_function{n,m}(Method.Iterations*Method.Deltat,-1i*FFTOperators1D.Gx); % Computing and storing the coupled time-dependent potential between components
    end
    % FOR each component where the stochastic potential is non null
    for m = FFTPhysics1D.StochasticPotential_function_Index{n}
        if (iscell(FFTPhysics1D.StochasticProcess_function))
            for m_noise = 1:length(FFTPhysics1D.StochasticProcess_function)
                FFTPhysics1D.StochasticProcess{m_noise} = (FFTPhysics1D.StochasticProcess_function{m_noise}(Method.Iterations*Method.Deltat,FFTGeometry1D.X)-FFTPhysics1D.StochasticProcess_function{m_noise}((Method.Iterations-1)*Method.Deltat,FFTGeometry1D.X))/Method.Deltat;
            end
        else
            FFTPhysics1D.StochasticProcess = (FFTPhysics1D.StochasticProcess_function(Method.Iterations*Method.Deltat,FFTGeometry1D.X) - FFTPhysics1D.StochasticProcess_function((Method.Iterations - 1)*Method.Deltat,FFTGeometry1D.X))/Method.Deltat;
        end
        FFTPhysics1D.StochasticPotential{n,m} = FFTPhysics1D.StochasticPotential_function{n,m}(FFTPhysics1D.StochasticProcess,FFTGeometry1D.X); % Computing and storing the stochastic potential
    end
    % FOR each component where the stochastic dispersion is non null
    for m = FFTPhysics1D.StochasticDispersion_function_Index{n}
        if (iscell(FFTPhysics1D.StochasticProcess_function))
            for m_noise = 1:length(FFTPhysics1D.StochasticProcess_function)
                FFTPhysics1D.StochasticProcess{m_noise} = (FFTPhysics1D.StochasticProcess_function{m_noise}(Method.Iterations*Method.Deltat,FFTGeometry1D.X)-FFTPhysics1D.StochasticProcess_function{m_noise}((Method.Iterations - 1)*Method.Deltat,FFTGeometry1D.X))/Method.Deltat;
            end
        else
            FFTPhysics1D.StochasticProcess = (FFTPhysics1D.StochasticProcess_function(Method.Iterations*Method.Deltat,FFTGeometry1D.X) - FFTPhysics1D.StochasticProcess_function((Method.Iterations - 1)*Method.Deltat,FFTGeometry1D.X))/Method.Deltat;
        end
        FFTPhysics1D.StochasticDispersion{n,m} = FFTPhysics1D.StochasticDispersion_function{n,m}(FFTPhysics1D.StochasticProcess,-1i*FFTOperators1D.Gx); % Computing and storing the stochastic potential
    end
end

%% Reshaping the ground states variables as a single vector for the iterative method
Phi_vect = zeros(Method.Ncomponents*FFTGeometry1D.Nx,1); % Initializing the vector to store each component's wave function
% FOR each component
for n = 1:Method.Ncomponents
    Phi_vect((1+(n-1)*FFTGeometry1D.Nx):(n*FFTGeometry1D.Nx)) = reshape(Phi{n},FFTGeometry1D.Nx,1); % Storing the component's wave function in the vector
end

%% Computing the wave functions after a single time step using Krylov iterative method for RSP
%If the computation is dynamic
if (strcmp(Method.Computation,'Dynamic'))
    Method.Deltat = 1i*Method.Deltat;
end
% IF one has chosen to use the Laplace preconditioner
if (strcmp(Method.Precond,'Laplace') == 1)
   % Using the iterative method which calls the operator with the 
   % Laplace preconditioner for the implicit part of the RSP scheme 
   % and the explicit part with the Laplace preconditioner applied. 
   % Moreover, it uses the initial components' wave functions as initial vector.
    [Phi_out, flag, relres, iter, resvec] = cgs(@(x)operator_Full_RSP_PLaplacian1d(x, Psi, FFTPsi, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D),...
    rhs_Full_RSP_PLaplacian1d(Phi_vect, Psi, FFTPsi, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D), Method.Iterative_tol, Method.Iterative_maxit,[],[],Phi_vect);
% ELSEIF one has chosen to use the Laplace preconditioner
elseif (strcmp(Method.Precond,'FLaplace') == 1)
   FFTPhysics1D.FPLaplace = RSPFPLaplace1d(Method, FFTPhysics1D, FFTGeometry1D);
   % Using the iterative method which calls the operator with the 
   % Laplace preconditioner for the implicit part of the RSP scheme 
   % and the explicit part with the Laplace preconditioner applied. 
   % Moreover, it uses the initial components' wave functions as initial vector.
    [Phi_out, flag, relres, iter, resvec] = cgs(@(x)operator_Full_RSP_FPLaplacian1d(x, Psi, FFTPsi, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D),...
    rhs_Full_RSP_FPLaplacian1d(Phi_vect, Psi, FFTPsi, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D), Method.Iterative_tol, Method.Iterative_maxit,[],[],Phi_vect);
elseif (strcmp(Method.Precond,'ThomasFermi') == 1)
   % Using the iterative method which calls the operator with the 
   % Thomas-Fermi preconditioner for the implicit part of the RSP scheme 
   % and the explicit part with the Thomas-Fermi preconditioner applied. 
   % Moreover, it uses the initial components' wave functions as initial vector.
    [Phi_out, flag, relres, iter, resvec] = cgs(@(x)operator_Full_RSP_PThomasFermi1d(x, Psi, FFTPsi, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D),...
    rhs_Full_RSP_PThomasFermi1d(Phi_vect, Psi, FFTPsi, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D), Method.Iterative_tol, Method.Iterative_maxit,[],[],Phi_vect);
elseif (strcmp(Method.Precond,'FThomasFermi') == 1)
   FFTPhysics1D.FPThomasFermi = RSPFPThomas_Fermi1d(Psi, FFTPsi, Method, FFTPhysics1D, FFTGeometry1D);
   % Using the iterative method which calls the operator with the full
   % Thomas-Fermi preconditioner for the implicit part of the RSP scheme 
   % and the explicit part with the Thomas-Fermi preconditioner applied. 
   % Moreover, it uses the initial components' wave functions as initial vector.
    [Phi_out, flag, relres, iter, resvec] = cgs(@(x)operator_Full_RSP_FPThomasFermi1d(x, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D),...
    rhs_Full_RSP_FPThomasFermi1d(Phi_vect, Psi, FFTPsi, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D), Method.Iterative_tol, Method.Iterative_maxit,[],[],Phi_vect);
elseif (strcmp(Method.Precond,'None') == 1)
   % Using the iterative method which calls the operator for the implicit
   % part of the CNSP scheme and the explicit part. 
   % Moreover, it uses the initial components' wave functions as initial vector.
   [Phi_out, flag, relres, iter, resvec] = cgs(@(x)operator_Full_RSP1d(x, Psi, FFTPsi, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D),...
    rhs_Full_RSP1d(Phi_vect, Psi, FFTPsi, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D), Method.Iterative_tol, Method.Iterative_maxit,[],[],Phi_vect);
end

%% Reshaping vector output from the iterative method as matrix
% FOR each component
for n = 1:Method.Ncomponents
    Phi{n} = reshape(Phi_out((1+(n-1)*FFTGeometry1D.Nx):(n*FFTGeometry1D.Nx)),FFTGeometry1D.Nx,1); % Storing the components' wave functions back in a cell array
end