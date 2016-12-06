%% Computation of a step in time using the full BESP-CNFG scheme
%% INPUTS:
%%          Phi: Initial wave functions in the 2D geometry for the FFT (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          FFTGeometry2D: Structure containing variables concerning the geometry of the problem in 2D in the FFT context (structure) (see FFTGeometry2D_Var2d.m)
%%          FFTPhysics2D: Structure containing variables concerning the physics of the problem in 2D in the FFT context (structure) (see FFTPhysics2D_Var2d.m)
%%          FFTOperators2D: Structure containing the derivative FFT operators (structure) (see FFTOperators2D_Var2d.m)
%% OUTPUTS:
%%          Phi: Wave functions computated with the BESP method on a single step (cell array)
%%          [flag,relres,iter,resvec]: Outputs from bicgstab (vector) (see bicgstab.m)
%% FUNCTIONS USED:
%%          bicgstab: To compute the wave functions by an iterative method (line 52,60 and 67)
%%          operator_Full_BESP_PLaplacian2d: To compute the implicit part of the BESP scheme with the Laplace preconditioner already applied (line 52)
%%          LinearLaplace_preconditioner2d: To apply the Laplace preconditioner to the explicit part of the BESP scheme (line 53)
%%          operator_Full_BESP_PThomasFermi2d: To compute the implicit part of the BESP scheme with the Thomas-Fermi preconditioner already applied (line 60)
%%          TF_preconditioner2d: To apply the Thomas-Fermi preconditioner to the explicit part of the BESP scheme (line 61)
%%          operator_Full_BESP2d: To compute the implicit part of the BESP scheme (line 67)

function [Phi, flag, relres, iter, resvec] = Local_Full_BESP_solution2d(Phi, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D)
%% Computing the explicit nonlinearity, non-local nonlinearity and time-dependent potential
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
    % FOR each component where the nonlinear gradient in the x direction is non null
    for m = FFTPhysics2D.GradientNLx_function_Index{n}
        FFTPhysics2D.GradientNLx{n,m} = FFTPhysics2D.GradientNLx_function{n,m}(Phi,FFTGeometry2D.X,FFTGeometry2D.Y,-1i*FFTOperators2D.Gx,-1i*FFTOperators2D.Gy); % Computing and storing the coupled nonlinear gradients in the x direction between components
    end
    % FOR each component where the nonlinear gradient in the y direction is non null
    for m = FFTPhysics2D.GradientNLy_function_Index{n}
        FFTPhysics2D.GradientNLy{n,m} = FFTPhysics2D.GradientNLy_function{n,m}(Phi,FFTGeometry2D.X,FFTGeometry2D.Y,-1i*FFTOperators2D.Gx,-1i*FFTOperators2D.Gy); % Computing and storing the coupled nonlinear gradients in the y direction between components
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
    
%% Computing the wave functions after a single time step using Krylov iterative method for BESP
%If the computation is dynamic
if (strcmp(Method.Computation,'Dynamic'))
    Method.Deltat = 1i*Method.Deltat;
end
% IF one has chosen to use the Laplace preconditioner
if (strcmp(Method.Precond,'Laplace') == 1)
   % Using the iterative method which calls the operator with the 
   % Laplace preconditionner for the implicit part of the BESP scheme 
   % and the explicit part with the Laplace preconditioner applied. 
   % Moreover, it uses the initial components' wave functions as initial vector.
   [Phi_out, flag, relres, iter, resvec] = bicgstab(@(x)operator_Full_BESP_PLaplacian2d(x, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D),...
   LinearLaplace_preconditioner2d((1/Method.Deltat)*Phi_vect, Method, FFTGeometry2D, FFTPhysics2D), Method.Iterative_tol, Method.Iterative_maxit,[],[],Phi_vect);
% ELSEIF one has chosen to use the full ThomasFermi preconditioner
elseif (strcmp(Method.Precond,'ThomasFermi') == 1)
   % Using the iterative method which calls the operator with the
   % Thomas-Fermi preconditioner for the implicit part of the BESP scheme 
   % and the explicit part with the Thomas-Fermi preconditioner applied. 
   % Moreover, it uses the initial components' wave functions as initial vector.
   [Phi_out, flag, relres, iter, resvec] = bicgstab(@(x)operator_Full_BESP_PThomasFermi2d(x, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D),...
   TF_preconditioner2d((1/Method.Deltat)*Phi_vect, Method, FFTGeometry2D, FFTPhysics2D),Method.Iterative_tol, Method.Iterative_maxit,[],[],Phi_vect);
% ELSEIF one has chosen to use the full ThomasFermi preconditioner
elseif (strcmp(Method.Precond,'FThomasFermi') == 1)
   FFTPhysics2D.FPThomasFermi = BESPFPThomas_Fermi2d(Method, FFTPhysics2D, FFTGeometry2D); % Computing the full Thomas Fermi preconditioner
   % Using the iterative method which calls the operator with the full
   % Thomas-Fermi preconditioner for the implicit part of the BESP scheme 
   % and the explicit part with the Thomas-Fermi preconditioner applied. 
   % Moreover, it uses the initial components' wave functions as initial vector.
   [Phi_out, flag, relres, iter, resvec] = bicgstab(@(x)operator_Full_BESP_FPThomasFermi2d(x, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D),...
   FTF_preconditioner2d((1/Method.Deltat)*Phi_vect, Method, FFTGeometry2D, FFTPhysics2D),Method.Iterative_tol, Method.Iterative_maxit,[],[],Phi_vect);
% ELSEIF one has chosen to use the full ThomasFermi preconditioner
elseif (strcmp(Method.Precond,'FLaplace') == 1)
   % Using the iterative method which calls the operator with the full
   % Laplace preconditioner for the implicit part of the BESP scheme 
   % and the explicit part with the Laplace preconditioner applied. 
   % Moreover, it uses the initial components' wave functions as initial vector.
   [Phi_out, flag, relres, iter, resvec] = bicgstab(@(x)operator_Full_BESP_FPLaplacian2d(x, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D),...
   FLP_preconditioner2d((1/Method.Deltat)*Phi_vect, Method, FFTGeometry2D, FFTPhysics2D),Method.Iterative_tol, Method.Iterative_maxit,[],[],Phi_vect);
% ELSEIF one has chosen not to use any preconditioner
elseif (strcmp(Method.Precond,'None') == 1)
   % Using the iterative method which calls the operator for the implicit 
   % part of the BESP scheme and the explicit part. Moreover, it uses the 
   % initial components' wave functions as initial vector.
   [Phi_out, flag, relres, iter, resvec] = bicgstab(@(x)operator_Full_BESP2d(x, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D),...
   (1/Method.Deltat)*Phi_vect, Method.Iterative_tol, Method.Iterative_maxit,[],[],Phi_vect);
end

%% Reshaping vector output from the iterative method as matrix
% FOR each component
for n = 1:Method.Ncomponents
    Phi{n} = reshape(Phi_out((1+(n-1)*FFTGeometry2D.N2):(n*FFTGeometry2D.N2)),FFTGeometry2D.Ny,FFTGeometry2D.Nx); % Storing the components' wave functions back in a cell array
end