%% Computation of a step in time using the RSP scheme
%% INPUT:
%%          Phi: Initial wave functions in the 3D geometry for the FFT (cell array)
%%          Psi: Relaxation variable corresponding to the nonlinearity (cell array)
%%          FFTPsi: Relaxation variable corresponding to the non local nonlinearity (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var3d.m)
%%          FFTGeometry3D: Structure containing variables concerning the geometry of the problem in 3D in the FFT context (structure) (see FFTGeometry3D_Var3d.m)
%%          FFTPhysics3D: Structure containing variables concerning the physics of the problem in 3D in the FFT context (structure) (see FFTPhysics3D_Var3d.m)
%%          FFTOperators3D: Structure containing the derivative FFT operators (structure) (see FFTOperators3D_Var3d.m)
%% OUTPUTS:
%%          Phi: Wave functions computated with the BESP method on a single step (cell array)
%%          [flag,relres,iter,resvec]: Outputs from bicgstab (vector) (see bicgstab.m)
%% FUNCTIONS USED:
%%          bicgstab: To compute the wave functions by an iterative method (line 50, 57 and 63)
%%          operator_Full_CNSP_PLaplacian3d: To compute the implicit part of the CNSP scheme with the Laplace preconditioner already applied (line 50)
%%          rhs_Full_CNSP_PLaplacian3d: To compute the explicit part of the CNSP scheme with the Laplace preconditioner already applied (line 51)
%%          operator_Full_CNSP_PThomasFermi3d: To compute the implicit part of the CNSP scheme with the Thomas-Fermi preconditioner already applied (line 57)
%%          rhs_Full_CNSP_PThomasFermi3d: To compute the explicit part of the CNSP scheme with the Thomas-Fermi preconditioner already applied (line 58)
%%          operator_Full_CNSP3d: To compute the implicit part of the CNSP scheme (line 63)
%%          rhs_Full_CNSP3d: To compute the explicit part of the CNSP scheme (line 64)

function [Phi, Psi, FFTPsi, flag, relres, iter, resvec] = Local_Full_RSP_solution3d(Phi, Psi, FFTPsi, Method, FFTGeometry3D, FFTPhysics3D, FFTOperators3D)

%% Computing the explicit nonlinearity, non local nonlinearity and the time-dependent potential
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component where the nonlinearity is non null
    for m = FFTPhysics3D.Nonlinearity_function_Index{n}
        FFTPhysics3D.Nonlinearity{n,m} = FFTPhysics3D.Nonlinearity_function{n,m}(Phi,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z); % Computing and storing the coupled nonlinearities between components
        Psi{n,m} = 2*FFTPhysics3D.Nonlinearity{n,m} - Psi{n,m}; % Computing the relaxation variable corresponding to the local nonlinearities
    end
    % FOR each component where the non-local nonlinearity is non null
    for m = FFTPhysics3D.FFTNonlinearity_function_Index{n}
        FFTPhysics3D.FFTNonlinearity{n,m} = FFTPhysics3D.FFTNonlinearity_function{n,m}(Phi,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z,-1i*FFTOperators3D.Gx,-1i*FFTOperators3D.Gy,-1i*FFTOperators3D.Gz); % Computing and storing the coupled non-local nonlinearities between components
        FFTPsi{n,m} = 2*FFTPhysics3D.FFTNonlinearity{n,m} - FFTPsi{n,m}; % Computing the relaxation variable corresponding to the local nonlinearities
    end
    % FOR each component where the time-dependent potential is non null
    for m = FFTPhysics3D.TimePotential_function_Index{n}
        FFTPhysics3D.TimePotentialImp{n,m} = FFTPhysics3D.TimePotential_function{n,m}(Method.Iterations*Method.Deltat,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z); % Computing and storing the implicit coupled time-dependent potentials between components
        FFTPhysics3D.TimePotentialExp{n,m} = FFTPhysics3D.TimePotential_function{n,m}((Method.Iterations - 1)*Method.Deltat,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z); % Computing and storing the explicit coupled time-dependent potentials between components
    end
     % FOR each component where the time-dependent dispersion is non null
    for m = FFTPhysics3D.TimeDispersion_function_Index{n}
        FFTPhysics3D.TimeDispersionImp{n,m} = FFTPhysics3D.TimeDispersion_function{n,m}(Method.Iterations*Method.Deltat,-1i*FFTOperators3D.Gx,-1i*FFTOperators3D.Gy,-1i*FFTOperators3D.Gz); % Computing and storing the implicit coupled time-dependent dispersion between components
        FFTPhysics3D.TimeDispersionExp{n,m} = FFTPhysics3D.TimeDispersion_function{n,m}((Method.Iterations - 1)*Method.Deltat,-1i*FFTOperators3D.Gx,-1i*FFTOperators3D.Gy,-1i*FFTOperators3D.Gz); % Computing and storing the explicit coupled time-dependent dispersion between components
    end
    % FOR each component where the stochastic potential is non null
    for m = FFTPhysics3D.StochasticPotential_function_Index{n}
        if (iscell(FFTPhysics3D.StochasticProcess_function))
            for m_noise = 1:length(FFTPhysics3D.StochasticProcess_function)
                FFTPhysics3D.StochasticProcess{m_noise} = (FFTPhysics3D.StochasticProcess_function{m_noise}(Method.Iterations*Method.Deltat,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z)-FFTPhysics3D.StochasticProcess_function{m_noise}((Method.Iterations - 1)*Method.Deltat,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z))/Method.Deltat;
            end
        else
            FFTPhysics3D.StochasticProcess = (FFTPhysics3D.StochasticProcess_function(Method.Iterations*Method.Deltat,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z)-FFTPhysics3D.StochasticProcess_function((Method.Iterations - 1)*Method.Deltat,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z))/Method.Deltat;
        end
        FFTPhysics3D.StochasticPotential{n,m} = FFTPhysics3D.StochasticPotential_function{n,m}(FFTPhysics3D.StochasticProcess,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z); % Computing and storing the stochastic potential
    end
    % FOR each component where the stochastic dispersion is non null
    for m = FFTPhysics3D.StochasticDispersion_function_Index{n}
        if (iscell(FFTPhysics3D.StochasticProcess_function))
            for m_noise = 1:length(FFTPhysics3D.StochasticProcess_function)
                FFTPhysics3D.StochasticProcess{m_noise} = (FFTPhysics3D.StochasticProcess_function{m_noise}(Method.Iterations*Method.Deltat,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z)-FFTPhysics3D.StochasticProcess_function{m_noise}((Method.Iterations - 1)*Method.Deltat,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z))/Method.Deltat;
            end
        else
            FFTPhysics3D.StochasticProcess = (FFTPhysics3D.StochasticProcess_function(Method.Iterations*Method.Deltat,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z)-FFTPhysics3D.StochasticProcess_function((Method.Iterations - 1)*Method.Deltat,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z))/Method.Deltat;
        end
        FFTPhysics3D.StochasticDispersion{n,m} = FFTPhysics3D.StochasticDispersion_function{n,m}(FFTPhysics3D.StochasticProcess,-1i*FFTOperators3D.Gx,-1i*FFTOperators3D.Gy,-1i*FFTOperators3D.Gz); % Computing and storing the stochastic potential
    end
end

%If the computation is dynamic
if (strcmp(Method.Computation,'Dynamic'))
    Method.Deltat = 1i*Method.Deltat;
end

%% Reshaping the ground states variables as a single vector for the iterative method
Phi_vect = zeros(Method.Ncomponents*FFTGeometry3D.N3,1); % Initializing the vector to store each component's wave function
% FOR each component
for n = 1:Method.Ncomponents
    Phi_vect((1+(n-1)*FFTGeometry3D.N3):(n*FFTGeometry3D.N3)) = reshape(Phi{n},FFTGeometry3D.N3,1); % Storing the component's wave function in the vector
end

%% Computing the wave functions after a single time step using Krylov iterative method for CNSP-CNFG
% IF one has chosen to use the Laplace preconditioner
if (strcmp(Method.Precond,'Laplace') == 1)
   % Using the iterative method which calls the operator with the 
   % Laplace preconditioner for the implicit part of the RSP scheme 
   % and the explicit part with the Laplace preconditioner applied. 
   % Moreover, it uses the initial components' wave functions as initial vector.
    [Phi_out, flag, relres, iter, resvec] = cgs(@(x)operator_Full_RSP_PLaplacian3d(x, Psi, FFTPsi, Method, FFTGeometry3D, FFTPhysics3D, FFTOperators3D),...
    rhs_Full_RSP_PLaplacian3d(Phi_vect, Psi, FFTPsi, Method, FFTGeometry3D, FFTPhysics3D, FFTOperators3D), Method.Iterative_tol, Method.Iterative_maxit,[],[],Phi_vect);
% IF one has chosen to use the Laplace preconditioner
elseif (strcmp(Method.Precond,'FLaplace') == 1)
   FFTPhysics3D.FPLaplace = RSPFPLaplace3d(Method, FFTPhysics3D, FFTGeometry3D);
   % Using the iterative method which calls the operator with the 
   % Laplace preconditioner for the implicit part of the RSP scheme 
   % and the explicit part with the Laplace preconditioner applied. 
   % Moreover, it uses the initial components' wave functions as initial vector.
    [Phi_out, flag, relres, iter, resvec] = cgs(@(x)operator_Full_RSP_FPLaplacian3d(x, Psi, FFTPsi, Method, FFTGeometry3D, FFTPhysics3D, FFTOperators3D),...
    rhs_Full_RSP_FPLaplacian3d(Phi_vect, Psi, FFTPsi, Method, FFTGeometry3D, FFTPhysics3D, FFTOperators3D), Method.Iterative_tol, Method.Iterative_maxit,[],[],Phi_vect);
% ELSEIF one has chosen to use the ThomasFermi preconditioner
elseif (strcmp(Method.Precond,'ThomasFermi') == 1)
   % Using the iterative method which calls the operator with the 
   % Thomas-Fermi preconditioner for the implicit part of the RSP scheme 
   % and the explicit part with the Thomas-Fermi preconditioner applied. 
   % Moreover, it uses the initial components' wave functions as initial vector.
    [Phi_out, flag, relres, iter, resvec] = cgs(@(x)operator_Full_RSP_PThomasFermi3d(x, Psi, FFTPsi, Method, FFTGeometry3D, FFTPhysics3D, FFTOperators3D),...
    rhs_Full_RSP_PThomasFermi3d(Phi_vect, Psi, FFTPsi, Method, FFTGeometry3D, FFTPhysics3D, FFTOperators3D), Method.Iterative_tol, Method.Iterative_maxit,[],[],Phi_vect);
% ELSEIF one has chosen to use the full ThomasFermi preconditioner
elseif (strcmp(Method.Precond,'FThomasFermi') == 1)
FFTPhysics3D.FPThomasFermi = RSPFPThomas_Fermi3d(Psi, FFTPsi, Method, FFTPhysics3D, FFTGeometry3D);
   % Using the iterative method which calls the operator with the 
   % Thomas-Fermi preconditioner for the implicit part of the RSP scheme 
   % and the explicit part with the Thomas-Fermi preconditioner applied. 
   % Moreover, it uses the initial components' wave functions as initial vector.
    [Phi_out, flag, relres, iter, resvec] = cgs(@(x)operator_Full_RSP_FPThomasFermi3d(x, Psi, FFTPsi, Method, FFTGeometry3D, FFTPhysics3D, FFTOperators3D),...
    rhs_Full_RSP_FPThomasFermi3d(Phi_vect, Psi, FFTPsi, Method, FFTGeometry3D, FFTPhysics3D, FFTOperators3D), Method.Iterative_tol, Method.Iterative_maxit,[],[],Phi_vect);
% ELSEIF one has chosen not to use any preconditioner
elseif (strcmp(Method.Precond,'None') == 1)
   % Using the iterative method which calls the operator for the implicit
   % part of the RSP scheme and the explicit part. 
   % Moreover, it uses the initial components' wave functions as initial vector.
   [Phi_out, flag, relres, iter, resvec] = cgs(@(x)operator_Full_RSP3d(x, Psi, FFTPsi, Method, FFTGeometry3D, FFTPhysics3D, FFTOperators3D),...
    rhs_Full_RSP3d(Phi_vect, Psi, FFTPsi, Method, FFTGeometry3D, FFTPhysics3D, FFTOperators3D), Method.Iterative_tol, Method.Iterative_maxit,[],[],Phi_vect);
end

%% Reshaping vector output from the iterative method as matrix
% FOR each component
for n = 1:Method.Ncomponents
    Phi{n} = reshape(Phi_out((1+(n-1)*FFTGeometry3D.N3):(n*FFTGeometry3D.N3)),FFTGeometry3D.Ny,FFTGeometry3D.Nx,FFTGeometry3D.Nz); % Storing the components' wave functions back in a cell array
end