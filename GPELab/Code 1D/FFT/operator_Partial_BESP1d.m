%% Solving BESP-CNFG scheme on a single step of time using a direct iterative method with Laplace preconditionner
%% INPUTS:
%%          Phi_in: Initial components' wave functions (cell array)
%%          Nonlinear_phi: Nonlinearities of each component (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var1d.m)
%%          FFTGeometry1D: Structure containing variables concerning the geometry of the problem in 1D in the FFT context (structure) (see FFTGeometry1D_Var1d.m)
%%          FFTPhysics1D: Structure containing variables concerning the physics of the problem in 1D in the FFT context (structure) (see FFTPhysics1D_Var1d.m)
%%          FFTOperators1D: Structure containing the derivative FFT operators (structure) (see FFTOperators1D_Var1d.m)
%% OUTPUTS:
%%          Phi_out: Components' wave functions with the operators applied (cell array)
%%          iter: Total number of iterations before convergence (double)

function [Phi_out,iter] = operator_Partial_BESP1d(Phi_in, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D)
%% Initialization of variables
Phi_tmp = Phi_in; % Storing the initial wave functions in a temporary variable
iter = 0; % Initialization of the iterations count
Phi_evo = 2*Method.Iterative_tol; % Initialization of the evolution criterion
Gradx_Phi = cell(1,Method.Ncomponents); % Initializing the variable for the gradient of the components' wave functions in the x direction
alpha = zeros(1,Method.Ncomponents); % Initializing the variable for the average between the maximum and the minimum of the potential and nonlinear terms
% FOR each component
for n = 1:Method.Ncomponents
b_min = min(FFTPhysics1D.Potential{n,n} + FFTPhysics1D.Beta*FFTPhysics1D.Nonlinearity{n,n}.*Phi_in{n}); % Computing the minimum of the potential and nonlinear terms
b_max = max(FFTPhysics1D.Potential{n,n} + FFTPhysics1D.Beta*FFTPhysics1D.Nonlinearity{n,n}.*Phi_in{n}); % Computing the maximum of the potential and nonlinear terms
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
    % IF there are gradient in the x direction functions for this
    % component to compute
    if (FFTPhysics1D.Gradientx_compute_Index{n})
        Gradx_Phi{n} = ifft(FFTOperators1D.Gx.*fft(Phi_tmp{n})); % First order derivating the component's wave function on the x direction via FFT and iFFT
    end
end
%% Applying the operators with the Laplace preconditionner
% FOR each component
for n = 1:Method.Ncomponents
    GPE_Phi = zeros(FFTGeometry1D.Nx,1); % Initializing the variable that will contain the wave function with the operators applied
    % FOR each component where the gradient in the x direction is non null
    for m = FFTPhysics1D.Gradientx_function_Index{n}
        GPE_Phi = GPE_Phi - FFTPhysics1D.Gradientx{n,m}.*Gradx_Phi{m}; % Computing the wave function with the gradient operators in the x direction applied
    end
    % FOR each component where the potential is non null
    for m = FFTPhysics1D.Potential_function_Index{n}
        GPE_Phi = GPE_Phi - FFTPhysics1D.Potential{n,m}.*Phi_tmp{m}; % Computing the wave function with the potential operator applied
    end
    % FOR each component where the nonlinearity is non null
    for m = FFTPhysics1D.Nonlinearity_function_Index{n}
        GPE_Phi = GPE_Phi - FFTPhysics1D.Beta*FFTPhysics1D.Nonlinearity{n,m}.*Phi_tmp{m}; % Computing the wave function with the nonlinearity operator applied
    end
    % FOR each component where the non-local nonlinearity is non null
    for m = FFTPhysics1D.Nonlinearity_function_Index{n}
        GPE_Phi = GPE_Phi - FFTPhysics1D.FFTNonlinearity{n,m}.*Phi_tmp{m}; % Computing the wave function with the non-local nonlinearity operator applied
    end
    % FOR each component where the time-dependent potential is non null
    for m = FFTPhysics1D.TimePotential_function_Index{n}
        GPE_Phi = GPE_Phi - FFTPhysics1D.TimePotential{n,m}.*Phi_tmp{m}; % Computing the wave function with the time-dependent potential operator applied
    end
    % FOR each component where the dispersion is non null
    for m = FFTPhysics1D.Dispersion_function_Index{n}
        % IF it is an extradiagonal term
        if (m ~= n)
        GPE_Phi = GPE_Phi - ifft(FFTPhysics1D.Dispersion{n,m}.*fft(Phi_tmp{m}));
        end
    end
        GPE_Phi = GPE_Phi + alpha(n).*Phi_tmp{n}; % Computing the wave function with the potential and nonlinear operators applied
        Phi_tmp{n} = ifft(1./(1/Method.Deltat+alpha(n) + FFTPhysics1D.Dispersion{n,n}).*fft(1/Method.Deltat*Phi_in{n} + GPE_Phi));  % Computing the IFFT of the explicit and implicit wave function then apply the Laplace preconditionner
        Phi_evo = max(Phi_evo,max(abs(Phi_tmp2{n} - Phi_tmp{n}))); % Computing the global evolution of the wave functions
end
    iter = iter +1; % Updating the number of iterations
end

%% Storing the output
Phi_out = Phi_tmp; % Storing the output