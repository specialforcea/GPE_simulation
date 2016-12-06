%% Applying the operators corresponding to the implicit part of the BESP-CNFG scheme with the Laplace preconditioner
%% INPUTS:
%%          Phi_in: Initial components' wave functions (vector)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          FFTGeometry2D: Structure containing variables concerning the geometry of the problem in 2D in the FFT context (structure) (see FFTGeometry2D_Var2d.m)
%%          FFTPhysics2D: Structure containing variables concerning the physics of the problem in 2D in the FFT context (structure) (see FFTPhysics2D_Var2d.m)
%%          FFTOperators2D: Structure containing the derivative FFT operators (structure) (see FFTOperators2D_Var2d.m)
%% OUTPUT:
%%          Phi_out: Components' wave functions with the operators applied (vector)

function [Phi_out] = operator_Full_BESP_PLaplacian2d(Phi_in, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D)
%% Initialization of variables
Phi = cell(1,Method.Ncomponents); % Initializing the variable for the components' wave functions
Gradx_Phi = cell(1,Method.Ncomponents); % Initializing the variable for the gradient of the components' wave functions in the x direction
Grady_Phi = cell(1,Method.Ncomponents); % Initializing the variable for the gradient of the components' wave functions in the y direction
Phi_out = zeros(FFTGeometry2D.Ny,Method.Ncomponents*FFTGeometry2D.Nx); % Initializing the variable for the components' wave functions with the operators and the preconditioner applied

%% Computing gradients
% FOR each component
for n = 1:Method.Ncomponents
    Phi{n} = reshape(Phi_in((1+(n-1)*FFTGeometry2D.N2):(n*FFTGeometry2D.N2)),FFTGeometry2D.Ny,FFTGeometry2D.Nx); % Storing the components' wave functions back in a cell array
    % IF there are gradient in the y direction functions for this
    % component to compute
    if (FFTPhysics2D.Gradienty_compute_Index{n})||(FFTPhysics2D.GradientNLy_compute_Index{n})
        Grady_Phi{n} = ifft(FFTOperators2D.Gy.*fft(Phi{n},[],1),[],1); % First order derivating the component's wave function on the y direction via FFT and iFFT
    end
    % IF there are gradient in the x direction functions for this
    % component to compute
    if (FFTPhysics2D.Gradientx_compute_Index{n})||(FFTPhysics2D.GradientNLx_compute_Index{n})
        Gradx_Phi{n} = ifft(FFTOperators2D.Gx.*fft(Phi{n},[],2),[],2); % First order derivating the component's wave function on the x direction via FFT and iFFT
    end
end

%% Applying the operators
% FOR each component
for n = 1:Method.Ncomponents
    GPE_Phi = zeros(FFTGeometry2D.Ny,FFTGeometry2D.Nx); % Initializing the variable that will contain the wave function with the operators applied
    % FOR each component where the gradient in the x direction is non null
    for m = FFTPhysics2D.Gradientx_function_Index{n}
        GPE_Phi = GPE_Phi + FFTPhysics2D.Gradientx{n,m}.*Gradx_Phi{m}; % Computing the wave function with the gradient operators in the x direction applied
    end
    % FOR each component where the gradient in the y direction is non null
    for m = FFTPhysics2D.Gradienty_function_Index{n}
        GPE_Phi = GPE_Phi + FFTPhysics2D.Gradienty{n,m}.*Grady_Phi{m}; % Computing the wave function with the gradient operators in the y direction applied
    end
    % FOR each component where the nonlinear gradient in the x direction is non null
    for m = FFTPhysics2D.GradientNLx_function_Index{n}
        GPE_Phi = GPE_Phi + FFTPhysics2D.GradientNLx{n,m}.*Gradx_Phi{m}; % Computing the wave function with the nonlinear gradient operators in the x direction applied
    end
    % FOR each component where the nonlinear gradient in the y direction is non null
    for m = FFTPhysics2D.GradientNLy_function_Index{n}
        GPE_Phi = GPE_Phi + FFTPhysics2D.GradientNLy{n,m}.*Grady_Phi{m}; % Computing the wave function with the nonlinear gradient operators in the y direction applied
    end
    % FOR each component where the potential is non null
    for m = FFTPhysics2D.Potential_function_Index{n}
        GPE_Phi = GPE_Phi + FFTPhysics2D.Potential{n,m}.*Phi{m}; % Computing the wave function with the potential operator applied
    end
    % FOR each component where the nonlinearity is non null
    for m = FFTPhysics2D.Nonlinearity_function_Index{n}
        GPE_Phi = GPE_Phi + FFTPhysics2D.Beta*FFTPhysics2D.Nonlinearity{n,m}.*Phi{m}; % Computing the wave function with the nonlinearity operator applied
    end
    % FOR each component where the non-local nonlinearity is non null
    for m = FFTPhysics2D.FFTNonlinearity_function_Index{n}
        GPE_Phi = GPE_Phi + FFTPhysics2D.Beta*FFTPhysics2D.FFTNonlinearity{n,m}.*Phi{m}; % Computing the wave function with the non-local nonlinearity operator applied
    end
    % FOR each component where the time-dependent potential is non null
    for m = FFTPhysics2D.TimePotential_function_Index{n}
        GPE_Phi = GPE_Phi + FFTPhysics2D.TimePotential{n,m}.*Phi{m}; % Computing the wave function with the time-dependent potential operator applied
    end
    % FOR each component where the dipersion is non null
    for m = FFTPhysics2D.Dispersion_function_Index{n}
        % IF it is an extradiagonal term
        if (m ~= n)
            GPE_Phi = GPE_Phi + ifft2(FFTPhysics2D.Dispersion{n,m}.*fft2(Phi{m})); % Computing the wave function with the dispersion operator applied
        end
    end
    GPE_Phi = Phi{n} + ifft2((1/Method.Deltat + FFTPhysics2D.Dispersion{n,n}).^(-1).*fft2(GPE_Phi)); % Computing the IFFT of the wave function with the operators applied then apply the Laplace preconditioner
    Phi_out(:,(1+(n-1)*FFTGeometry2D.Nx):(n*FFTGeometry2D.Nx)) = GPE_Phi; % Storing the wave function with the operators and the preconditionner applied
end

%% Reshapping as a vector the output
Phi_out = reshape(Phi_out,Method.Ncomponents*FFTGeometry2D.N2,1); % Reshapping the wave functions as a vector