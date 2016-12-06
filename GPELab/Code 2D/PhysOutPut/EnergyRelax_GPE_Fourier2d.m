%% Computation of the relaxation energy of the GPE using the FFT operators
%% INPUTS:
%%          Phi: Wave functions (cell array)
%%          Psi: Relaxation functions (cell array)
%%          Psitmp: Temporary relaxation functions (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          FFTGeometry2D: Structure containing variables concerning the geometry of the problem in 2D in the FFT context (structure) (see FFTGeometry2D_Var2d.m)
%%          FFTPhysics2D: Structure containing variables concerning the physics of the problem in 2D in the FFT context (structure) (see FFTPhysics2D_Var2d.m)
%%          FFTOperators2D: Structure containing the derivative FFT operators (structure) (see FFTOperators2D_Var2d.m)
%% OUTPUT:
%%          Energy: Energy of the wave functions (cell array)
%% FUNCTIONS USED:
%%          L2_norm2d: To integrate the local energy (line 62)

function [Energy] = EnergyRelax_GPE_Fourier2d(Phi, Psi, Psitmp, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D)
%% Initialization
Gradx_Phi = cell(1,Method.Ncomponents); % Initializing the variable for the gradient of the components' wave functions in the x direction
Grady_Phi = cell(1,Method.Ncomponents); % Initializing the variable for the gradient of the components' wave functions in the y direction

%% Computation of the gradients
% FOR each component
for n = 1:Method.Ncomponents
        Grady_Phi{n} = ifft(FFTOperators2D.Gy.*fft(Phi{n},[],1),[],1); % First order derivating the component's wave function on the y direction via FFT and iFFT
        Gradx_Phi{n} = ifft(FFTOperators2D.Gx.*fft(Phi{n},[],2),[],2); % First order derivating the component's wave function on the x direction via FFT and iFFT
end

%% Computation of the energy
% FOR each component
for n = 1:Method.Ncomponents
    %% Initialization 
    Transport_energy = zeros(FFTGeometry2D.Ny,FFTGeometry2D.Nx); % Initialization of the local transport energy
    Nonlinear_energy = zeros(FFTGeometry2D.Ny,FFTGeometry2D.Nx); % Initialization of the local nonlinear energy
    Potential_energy = zeros(FFTGeometry2D.Ny,FFTGeometry2D.Nx); % Initialization of the local potential energy
    %% Computation of the transport, nonlinear and potential local energies
    % FOR each component where the gradient in the x direction is non null
    for m = FFTPhysics2D.Gradientx_function_Index{n}
        Transport_energy = Transport_energy + 0*real(FFTPhysics2D.Gradientx_function{n,m}(FFTGeometry2D.X,FFTGeometry2D.Y).*Gradx_Phi{m}.*conj(Phi{n})); % Computing the local potential energy for the gradient in the x direction for the component
    end
    % FOR each component where the gradient in the y direction is non null
    for m = FFTPhysics2D.Gradienty_function_Index{n}
        Transport_energy = Transport_energy + 0*real(FFTPhysics2D.Gradienty_function{n,m}(FFTGeometry2D.X,FFTGeometry2D.Y).*Grady_Phi{m}.*conj(Phi{n})); % Computing the local potential energy for the gradient in the y direction for the component
    end
    % FOR each component where the potential is non null
    for m = FFTPhysics2D.Potential_function_Index{n}
        Potential_energy = Potential_energy + real(FFTPhysics2D.Potential_function{n,m}(FFTGeometry2D.X,FFTGeometry2D.Y).*Phi{m}.*conj(Phi{n})); % Computing the local potential energy for the component
    end
    % FOR each component where the time-dependent potential is non null
    for m = FFTPhysics2D.TimePotential_function_Index{n}
        Potential_energy = Potential_energy + real(FFTPhysics2D.TimePotential_function{n,m}(Method.Iterations*Method.Deltat,FFTGeometry2D.X,FFTGeometry2D.Y).*Phi{m}.*conj(Phi{n})); % Computing and storing the local time-dependent potential energy for the component
    end
    % FOR each component where the nonlinearity is non null
    for m = FFTPhysics2D.Nonlinearity_function_Index{n}
        % IF the energy's nonlinearity is defined
        if (isempty(FFTPhysics2D.Nonlinearity_energy_function{n,m}) == 0)
            Nonlinear_energy = Nonlinear_energy + real(FFTPhysics2D.Beta*FFTPhysics2D.Nonlinearity_energy_function{n,m}(Phi,FFTGeometry2D.X,FFTGeometry2D.Y).*Phi{m}.*conj(Phi{n})); % Computing of the local nonlinear energy
        % ELSE if the energy's nonlinearity is not defined
        else
            Nonlinear_energy = Nonlinear_energy + (FFTPhysics2D.Beta/2)*real(Psi{n,m}.*Psitmp{n,m}); % Computing of the local nonlinear energy
        end
    end
    % FOR each component where the nonlinearity is non null
    for m = FFTPhysics2D.FFTNonlinearity_function_Index{n}
        % IF the energy's nonlinearity is defined
        if (isempty(FFTPhysics2D.FFTNonlinearity_energy_function{n,m}) == 0)
            Nonlinear_energy = Nonlinear_energy + real(FFTPhysics2D.FFTNonlinearity_energy_function{n,m}(Phi,FFTGeometry2D.X,FFTGeometry2D.Y,-1i*FFTOperators2D.Gx,-1i*FFTOperators2D.Gy).*Phi{m}.*conj(Phi{n})); % Computing of the local nonlinear energy
        % ELSE if the energy's nonlinearity is not defined
        else
            Nonlinear_energy = Nonlinear_energy + real(FFTPhysics2D.FFTNonlinearity_function{n,m}(Phi,FFTGeometry2D.X,FFTGeometry2D.Y,-1i*FFTOperators2D.Gx,-1i*FFTOperators2D.Gy).*Phi{m}.*conj(Phi{n})); % Computing of the local nonlinear energy
        end
    end
    %% Computation of the local energy
    Local_energy{n} = FFTPhysics2D.Delta*(abs(Gradx_Phi{n}).^2 +abs(Grady_Phi{n}).^2) + Potential_energy + Nonlinear_energy + Transport_energy ; % Computing local energy
    %% Integration over the space
    Energy{n} = L2_norm2d(sqrt(Local_energy{n}), FFTGeometry2D)^2; % Integration of the local energy
end