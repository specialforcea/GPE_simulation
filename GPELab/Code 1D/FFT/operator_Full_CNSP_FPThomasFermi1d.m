%% Applying the operators corresponding to the implicit part of the CNSP-CNFG scheme with the full Thomas-Fermi preconditionner
%% INPUTS:
%%          Phi_in: Initial components' wave functions (vector)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var1d.m)
%%          FFTGeometry1D: Structure containing variables concerning the geometry of the problem in 1D in the FFT context (structure) (see FFTGeometry1D_Var1d.m)
%%          FFTPhysics1D: Structure containing variables concerning the physics of the problem in 1D in the FFT context (structure) (see FFTPhysics1D_Var1d.m)
%%          FFTOperators1D: Structure containing the derivative FFT operators (structure) (see FFTOperators1D_Var1d.m)
%% OUTPUT:
%%          Phi_out: Components' wave functions with the operators applied (vector)

function [Phi_out] = operator_Full_CNSP_FPThomasFermi1d(Phi_in, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D)
%% Initialization of variables
Phi = cell(1,Method.Ncomponents); % Initializing the variable for the components' wave functions
Gradx_Phi = cell(1,Method.Ncomponents); % Initializing the variable for the gradient of the components' wave functions in the x direction
Phi_out = zeros(Method.Ncomponents*FFTGeometry1D.Nx,1); % Initializing the variable for the components' wave functions with the operators and the preconditioner applied

%% Computing gradients
% FOR each component
for n = 1:Method.Ncomponents
    Phi{n} = reshape(Phi_in((1+(n-1)*FFTGeometry1D.Nx):(n*FFTGeometry1D.Nx)),FFTGeometry1D.Nx,1); % Storing the components' wave functions back in a cell array
    % IF there are gradient in the x direction functions for this
    % component to compute
    if (FFTPhysics1D.Gradientx_compute_Index{n})
        Gradx_Phi{n} = ifft(FFTOperators1D.Gx.*fft(Phi{n})); % First order derivating the component's wave function on the x direction via FFT and iFFT
    end
end

%% Applying the operators
% FOR each component
for n = 1:Method.Ncomponents
    GPE_Phi{n} = zeros(FFTGeometry1D.Nx,1); % Initializing the variable that will contain the wave function with the operators applied
    % FOR each component where the gradient in the x direction is non null
    for m = FFTPhysics1D.Gradientx_function_Index{n}
        GPE_Phi{n} = GPE_Phi{n} + FFTPhysics1D.Gradientx{n,m}.*Gradx_Phi{m}/2; % Computing the wave function with the gradient operators in the x direction applied
    end
    % FOR each component where the dispersion is non null
    for m = FFTPhysics1D.Dispersion_function_Index{n}
        % IF it is an extradiagonal term
        if (m ~= n)
        GPE_Phi{n} = GPE_Phi{n} + ifft(FFTPhysics1D.Dispersion{n,m}.*fft(Phi{m}))/2;
        end
    end
    GPE_Phi{n} = GPE_Phi{n} + ifft(FFTPhysics1D.Dispersion{n,n}/2.*fft(Phi{n})); % Computing the wave function with the Laplacian operator applied
end

%% Applying the full Thomas Fermi preconditioner
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
    Phi{n} = Phi{n} + FFTPhysics1D.FPThomasFermi{n,m}.*GPE_Phi{m}; % Applying the Thomas-Fermi preconditioner
    end
    Phi_out((1+(n-1)*FFTGeometry1D.Nx):(n*FFTGeometry1D.Nx)) = Phi{n}; % Storing the wave function with the operators and the preconditioner applied
end

%% Reshapping as a vector the output
Phi_out = reshape(Phi_out,Method.Ncomponents*FFTGeometry1D.Nx,1); % Reshapping the wave functions as a vector