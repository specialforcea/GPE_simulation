%% Computation of a step in time using the full BESP scheme
%% INPUTS:
%%          Phi: Initial wave functions in the 1D geometry for the FFT (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var1d.m)
%%          FFTGeometry1D: Structure containing variables concerning the geometry of the problem in 1D in the FFT context (structure) (see FFTGeometry1D_Var1d.m)
%%          FFTPhysics1D: Structure containing variables concerning the physics of the problem in 1D in the FFT context (structure) (see FFTPhysics1D_Var1d.m)
%%          FFTOperators1D: Structure containing the derivative FFT operators (structure) (see FFTOperators1D_Var1d.m)
%% OUTPUTS:
%%          Phi: Wave functions computated with the BESP method on a single step (cell array)
%%          [flag,relres,iter,resvec]: Outputs from bicgstab (vector) (see bicgstab.m)
%% FUNCTIONS USED:
%%          bicgstab: To compute the wave functions by an iterative method (line 52,60 and 67)
%%          operator_Full_BESP_PLaplacian1d: To compute the implicit part of the BESP scheme with the Laplace preconditioner already applied (line 52)
%%          LinearLaplace_preconditioner1d: To apply the Laplace preconditioner to the explicit part of the BESP scheme (line 53)
%%          operator_Full_BESP_PThomasFermi1d: To compute the implicit part of the BESP scheme with the Thomas-Fermi preconditioner already applied (line 60)
%%          TF_preconditioner1d: To apply the Thomas-Fermi preconditioner to the explicit part of the BESP scheme (line 61)
%%          operator_Full_BESP1d: To compute the implicit part of the BESP scheme (line 67)

function [Phi, flag, relres, iter, resvec] = Local_Full_BESP_solution1d(Phi, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D)

%% Computing the explicit nonlinearity, the explicit non-local nonlinearity and the time-dependent potential
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

%% Reshaping the ground states variables as a single vector for the iterative method
Phi_vect = zeros(Method.Ncomponents*FFTGeometry1D.Nx,1); % Initializing the vector to store each component's wave function
% FOR each component
for n = 1:Method.Ncomponents
    Phi_vect((1+(n-1)*FFTGeometry1D.Nx):(n*FFTGeometry1D.Nx)) = reshape(Phi{n},FFTGeometry1D.Nx,1); % Storing the component's wave function in the vector
end

%% Computing the wave functions after a single time step using Krylov iterative method for BESP
%If the computation is dynamic
if (strcmp(Method.Computation,'Dynamic'))
    Method.Deltat = 1i*Method.Deltat;
end
% IF one has chosen to use the Laplace preconditioner
if (strcmp(Method.Precond,'Laplace') == 1)
   % Using the iterative method which calls the operator with the 
   % Laplace preconditioner for the implicit part of the BESP scheme 
   % and the explicit part with the Laplace preconditioner applied. 
   % Moreover, it uses the initial components' wave functions as initial vector.
   [Phi_out, flag, relres, iter, resvec] = bicgstab(@(x)operator_Full_BESP_PLaplacian1d(x, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D),...
   LinearLaplace_preconditioner1d((1/Method.Deltat)*Phi_vect, Method, FFTGeometry1D, FFTPhysics1D), Method.Iterative_tol, Method.Iterative_maxit,[],[],Phi_vect);
% ELSEIF one has chosen to use the ThomasFermi preconditioner
elseif (strcmp(Method.Precond,'ThomasFermi') == 1)
   % Using the iterative method which calls the operator with the 
   % Thomas-Fermi preconditioner for the implicit part of the BESP scheme 
   % and the explicit part with the Thomas-Fermi preconditioner applied. 
   % Moreover, it uses the initial components' wave functions as initial vector.
   [Phi_out, flag, relres, iter, resvec] = bicgstab(@(x)operator_Full_BESP_PThomasFermi1d(x, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D),...
   TF_preconditioner1d((1/Method.Deltat)*Phi_vect, Method, FFTGeometry1D, FFTPhysics1D),Method.Iterative_tol, Method.Iterative_maxit,[],[],Phi_vect);
% ELSEIF one has chosen to use the full ThomasFermi preconditioner
elseif (strcmp(Method.Precond,'FThomasFermi') == 1)
    FFTPhysics1D.FPThomasFermi = BESPFPThomas_Fermi1d(Method, FFTPhysics1D, FFTGeometry1D);
   % Using the iterative method which calls the operator with the full
   % Thomas-Fermi preconditioner for the implicit part of the BESP scheme 
   % and the explicit part with the Thomas-Fermi preconditioner applied. 
   % Moreover, it uses the initial components' wave functions as initial vector.
   [Phi_out, flag, relres, iter, resvec] = bicgstab(@(x)operator_Full_BESP_FPThomasFermi1d(x, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D),...
   FTF_preconditioner1d((1/Method.Deltat)*Phi_vect, Method, FFTGeometry1D, FFTPhysics1D),Method.Iterative_tol, Method.Iterative_maxit,[],[],Phi_vect);
% ELSEIF one has chosen to use the full Laplace preconditioner
elseif (strcmp(Method.Precond,'FLaplace') == 1)
   % Using the iterative method which calls the operator with the full
   % Thomas-Fermi preconditioner for the implicit part of the BESP scheme 
   % and the explicit part with the Laplace preconditioner applied. 
   % Moreover, it uses the initial components' wave functions as initial vector.
   [Phi_out, flag, relres, iter, resvec] = bicgstab(@(x)operator_Full_BESP_FPLaplacian1d(x, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D),...
   FLP_preconditioner1d((1/Method.Deltat)*Phi_vect, Method, FFTGeometry1D, FFTPhysics1D),Method.Iterative_tol, Method.Iterative_maxit,[],[],Phi_vect);
% ELSEIF one has chosen not to use any preconditioner
elseif (strcmp(Method.Precond,'None') == 1)
   % Using the iterative method which calls the operator for the implicit 
   % part of the BESP scheme and the explicit part. Moreover, it uses the 
   % initial components' wave functions as initial vector.
   [Phi_out, flag, relres, iter, resvec] = bicgstab(@(x)operator_Full_BESP1d(x, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D),...
   (1/Method.Deltat)*Phi_vect, Method.Iterative_tol, Method.Iterative_maxit,[],[],Phi_vect);
end

%% Reshaping vector output from the iterative method as matrix
% FOR each component
for n = 1:Method.Ncomponents
    Phi{n} = reshape(Phi_out((1+(n-1)*FFTGeometry1D.Nx):(n*FFTGeometry1D.Nx)),FFTGeometry1D.Nx,1); % Storing the components' wave functions back in a cell array
end