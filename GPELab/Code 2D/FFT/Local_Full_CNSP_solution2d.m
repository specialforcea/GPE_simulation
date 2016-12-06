%% Computation of a step in time using the CNSP scheme
%% INPUT:
%%          Phi: Initial wave functions in the 2D geometry for the FFT (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          FFTGeometry2D: Structure containing variables concerning the geometry of the problem in 2D in the FFT context (structure) (see FFTGeometry2D_Var2d.m)
%%          FFTPhysics2D: Structure containing variables concerning the physics of the problem in 2D in the FFT context (structure) (see FFTPhysics2D_Var2d.m)
%%          FFTOperators2D: Structure containing the derivative FFT operators (structure) (see FFTOperators2D_Var2d.m)
%% OUTPUTS:
%%          Phi: Wave functions computated with the BESP method on a single step (cell array)
%%          [flag,relres,iter,resvec]: Outputs from bicgstab (vector) (see bicgstab.m)
%% FUNCTIONS USED:
%%          bicgstab: To compute the wave functions by an iterative method (line 53, 60 and 66)
%%          operator_Full_CNSP_PLaplacian2d: To compute the implicit part of the CNSP scheme with the Laplace preconditioner already applied (line 53)
%%          rhs_Full_CNSP_PLaplacian2d: To compute the explicit part of the CNSP scheme with the Laplace preconditioner already applied (line 54)
%%          operator_Full_CNSP_PThomasFermi2d: To compute the implicit part of the CNSP scheme with the Thomas-Fermi preconditioner already applied (line 60)
%%          rhs_Full_CNSP_PThomasFermi2d: To compute the explicit part of the CNSP scheme with the Thomas-Fermi preconditioner already applied (line 61)
%%          operator_Full_CNSP2d: To compute the implicit part of the CNSP scheme (line 66)
%%          rhs_Full_CNSP2d: To compute the explicit part of the CNSP scheme (line 67)

function [Phi,flag, relres, iter, resvec] = Local_Full_CNSP_solution2d(Phi, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D)

%% Computing the explicit nonlinearity, non-local nonlinearity and the time-depend potential
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

%% Reshaping the ground states variables as a single vector for the iterative method
Phi_vect = zeros(Method.Ncomponents*FFTGeometry2D.N2,1); % Initializing the vector to store each component's wave function
% FOR each component
for n = 1:Method.Ncomponents
    Phi_vect((1+(n-1)*FFTGeometry2D.N2):(n*FFTGeometry2D.N2)) = reshape(Phi{n},FFTGeometry2D.N2,1); % Storing the component's wave function in the vector
end

%% Computing the wave functions after a single time step using Krylov iterative method for CNSP
%If the computation is dynamic
if (strcmp(Method.Computation,'Dynamic'))
    Method.Deltat = 1i*Method.Deltat;
end
% IF one has chosen to use the Laplace preconditioner
if (strcmp(Method.Precond,'Laplace') == 1)
   % Using the iterative method which calls the operator with the 
   % Laplace preconditioner for the implicit part of the CNSP scheme 
   % and the explicit part with the Laplace preconditioner applied. 
   % Moreover, it uses the initial components' wave functions as initial vector.
    [Phi_out, flag, relres, iter, resvec] = bicgstab(@(x)operator_Full_CNSP_PLaplacian2d(x, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D),...
    rhs_Full_CNSP_PLaplacian2d(Phi_vect, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D), Method.Iterative_tol, Method.Iterative_maxit,[],[],Phi_vect);
elseif (strcmp(Method.Precond,'ThomasFermi') == 1)
   % Using the iterative method which calls the operator with the 
   % Thomas-Fermi preconditioner for the implicit part of the CNSP scheme 
   % and the explicit part with the Thomas-Fermi preconditioner applied. 
   % Moreover, it uses the initial components' wave functions as initial vector.
    [Phi_out, flag, relres, iter, resvec] = bicgstab(@(x)operator_Full_CNSP_PThomasFermi2d(x, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D),...
    rhs_Full_CNSP_PThomasFermi2d(Phi_vect, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D), Method.Iterative_tol, Method.Iterative_maxit,[],[],Phi_vect);
elseif (strcmp(Method.Precond,'FLaplace') == 1)
   % Using the iterative method which calls the operator with the 
   % Thomas-Fermi preconditioner for the implicit part of the CNSP scheme 
   % and the explicit part with the full Laplace preconditioner applied. 
   % Moreover, it uses the initial components' wave functions as initial vector.
    [Phi_out, flag, relres, iter, resvec] = bicgstab(@(x)operator_Full_CNSP_FPLaplacian2d(x, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D),...
    rhs_Full_CNSP_FPLaplacian2d(Phi_vect, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D), Method.Iterative_tol, Method.Iterative_maxit,[],[],Phi_vect);
elseif (strcmp(Method.Precond,'None') == 1)
   % Using the iterative method which calls the operator for the implicit
   % part of the CNSP scheme and the explicit part. 
   % Moreover, it uses the initial components' wave functions as initial vector.
   [Phi_out, flag, relres, iter, resvec] = bicgstab(@(x)operator_Full_CNSP2d(x, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D),...
    rhs_Full_CNSP2d(Phi_vect, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D), Method.Iterative_tol, Method.Iterative_maxit,[],[],Phi_vect);
end

%% Reshaping vector output from the iterative method as matrix
% FOR each component
for n = 1:Method.Ncomponents
    Phi{n} = reshape(Phi_out((1+(n-1)*FFTGeometry2D.N2):(n*FFTGeometry2D.N2)),FFTGeometry2D.Ny,FFTGeometry2D.Nx); % Storing the components' wave functions back in a cell array
end