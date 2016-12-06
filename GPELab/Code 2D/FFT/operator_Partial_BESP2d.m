%% Solving BESP-CNFG scheme on a single step of time using a direct iterative method with Laplace preconditionner
%% INPUTS:
%%          Phi_in: Initial components' wave functions (cell array)
%%          Nonlinear_phi: Nonlinearities of each component (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          FFTGeometry2D: Structure containing variables concerning the geometry of the problem in 2D in the FFT context (structure) (see FFTGeometry2D_Var2d.m)
%%          FFTPhysics2D: Structure containing variables concerning the physics of the problem in 2D in the FFT context (structure) (see FFTPhysics2D_Var2d.m)
%%          FFTOperators2D: Structure containing the derivative FFT operators (structure) (see FFTOperators2D_Var2d.m)
%% OUTPUTS:
%%          Phi_out: Components' wave functions with the operators applied (cell array)
%%          iter: Total number of iterations before convergence (double)

function [Phi_out,iter] = operator_Partial_BESP2d(Phi_in, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D)
%% Initialization of variables
Phi_tmp = Phi_in; % Storing the initial wave functions in a temporary variable
iter = 0; % Initialization of the iterations count
Phi_evo = 2*Method.Iterative_tol; % Initialization of the evolution criterion
Gradx_Phi = cell(1,Method.Ncomponents); % Initializing the variable for the gradient of the components' wave functions in the x direction
Grady_Phi = cell(1,Method.Ncomponents); % Initializing the variable for the gradient of the components' wave functions in the y direction
alpha = zeros(1,Method.Ncomponents); % Initializing the variable for the average between the maximum and the minimum of the potential and nonlinear terms
% FOR each component
for n = 1:Method.Ncomponents
b_min = min(min(FFTPhysics2D.Potential{n,n} + FFTPhysics2D.Beta*FFTPhysics2D.Nonlinearity{n,n}.*Phi_in{n})); % Computing the minimum of the potential and nonlinear terms
b_max = max(max(FFTPhysics2D.Potential{n,n} + FFTPhysics2D.Beta*FFTPhysics2D.Nonlinearity{n,n}.*Phi_in{n})); % Computing the maximum of the potential and nonlinear terms
alpha(n) = (b_max+b_min)/2; % Average between the maximum and the minimum of the potential and nonlinear terms
end

%% Computing the wave functions after a single time step using the direct iterative method for BESP-CNFG
% Stopping criterions:- the global evolution of the wave functions (maximum
%                     of each local evolution) must be lower than the 
%                     stopping criterion Method.Iterative_tol
%                     - the number of iterations must not exceed the
%                     maximum number of iterations Method.Iterative_maxit
while (Phi_evo >Method.Iterative_tol)&&(iter<Method.Iterative_maxit)
    Phi_evo = 0; % Reinitialization of the evolution criterion
    Phi_tmp2 = Phi_tmp; % Storing the temporary wave functions to compute the evolution criterion
%% Computing gradients
% FOR each component
for n = 1:Method.Ncomponents
    % IF there are gradient in the y direction functions for this
    % component to compute
    if (FFTPhysics2D.Gradienty_compute_Index{n})
        Grady_Phi{n} = ifft(FFTOperators2D.Gy.*fft(Phi_tmp{n},[],1),[],1); % First order derivating the component's wave function on the y direction via FFT and iFFT
    end
    % IF there are gradient in the x direction functions for this
    % component to compute
    if (FFTPhysics2D.Gradientx_compute_Index{n})
        Gradx_Phi{n} = ifft(FFTOperators2D.Gx.*fft(Phi_tmp{n},[],2),[],2); % First order derivating the component's wave function on the x direction via FFT and iFFT
    end
end
%% Applying the operators with the Laplace preconditionner
% FOR each component
for n = 1:Method.Ncomponents
    GPE_Phi = zeros(FFTGeometry2D.Ny,FFTGeometry2D.Nx); % Initializing the variable that will contain the wave function with the operators applied
    % FOR each component where the gradient in the x direction is non null
    for m = FFTPhysics2D.Gradientx_function_Index{n}
        GPE_Phi = GPE_Phi - FFTPhysics2D.Gradientx{n,m}.*Gradx_Phi{m}; % Computing the wave function with the gradient operators in the x direction applied
    end
    % FOR each component where the gradient in the y direction is non null
    for m = FFTPhysics2D.Gradienty_function_Index{n}
        GPE_Phi = GPE_Phi - FFTPhysics2D.Gradienty{n,m}.*Grady_Phi{m}; % Computing the wave function with the gradient operators in the y direction applied
    end
    % FOR each component where the potential is non null
    for m = FFTPhysics2D.Potential_function_Index{n}
        GPE_Phi = GPE_Phi - FFTPhysics2D.Potential{n,m}.*Phi_tmp{m}; % Computing the wave function with the potential operator applied
    end
    % FOR each component where the nonlinearity is non null
    for m = FFTPhysics2D.Nonlinearity_function_Index{n}
        GPE_Phi = GPE_Phi - FFTPhysics2D.Beta*FFTPhysics2D.Nonlinearity{n,m}.*Phi_tmp{m}; % Computing the wave function with the nonlinearity operator applied
    end
     % FOR each component where the nonlinearity is non null
    for m = FFTPhysics2D.FFTNonlinearity_function_Index{n}
        GPE_Phi = GPE_Phi - FFTPhysics2D.FFTNonlinearity{n,m}.*Phi_tmp{m}; % Computing the wave function with the nonlinearity operator applied
    end
    % FOR each component where the time-dependent potential is non null
    for m = FFTPhysics2D.TimePotential_function_Index{n}
        GPE_Phi = GPE_Phi - FFTPhysics2D.TimePotential{n,m}.*Phi_tmp{m}; % Computing the wave function with the time-dependent potential operator applied
    end
    % FOR each component where the dipersion is non null
    for m = FFTPhysics2D.Dispersion_function_Index{n}
        % IF it is an extradiagonal term
        if (m ~= n)
            GPE_Phi = GPE_Phi - ifft2(FFTPhysics2D.Dispersion{n,m}.*fft2(Phi_tmp{m})); % Computing the wave function with the dispersion operator applied
        end
    end
        GPE_Phi = GPE_Phi + alpha(n).*Phi_tmp{n}; % Computing the wave function with the potential and nonlinear operators applied
        Phi_tmp{n} = ifft2(1./(1/Method.Deltat+alpha(n)-FFTPhysics2D.Delta*(FFTOperators2D.Dx+FFTOperators2D.Dy) + FFTPhysics2D.Dispersion{n,n}).*fft2(1/Method.Deltat*Phi_in{n} + GPE_Phi));  % Computing the IFFT of the explicit and implicit wave function then apply the Laplace preconditionner
        Phi_evo = max(Phi_evo,max(max(abs(Phi_tmp2{n} - Phi_tmp{n})))); % Computing the global evolution of the wave functions
end
    iter = iter +1; % Updating the number of iterations
end

%% Storing the output
Phi_out = Phi_tmp; % Storing the output