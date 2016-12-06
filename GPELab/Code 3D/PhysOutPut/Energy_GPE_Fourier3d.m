%% Computation of the energy of the GPE using the FFT operators
%% INPUTS:
%%          Phi: Wave functions (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var3d.m)
%%          FFTGeometry3D: Structure containing variables concerning the geometry of the problem in 3D in the FFT context (structure) (see FFTGeometry3D_Var3d.m)
%%          FFTPhysics3D: Structure containing variables concerning the physics of the problem in 3D in the FFT context (structure) (see FFTPhysics3D_Var3d.m)
%%          FFTOperators3D: Structure containing the derivative FFT operators (structure) (see FFTOperators3D_Var3d.m)
%% OUTPUT:
%%          Energy: Energy of the wave functions (cell array)
%% FUNCTIONS USED:
%%          L2_norm3d: To integrate the local energy (line 62)

function [Energy] = Energy_GPE_Fourier3d(Phi, Method, FFTGeometry3D, FFTPhysics3D, FFTOperators3D)
%% Initialization
Gradx_Phi = cell(1,Method.Ncomponents); % Initializing the variable for the gradient of the components' wave functions in the x direction
Grady_Phi = cell(1,Method.Ncomponents); % Initializing the variable for the gradient of the components' wave functions in the y direction
Gradz_Phi = cell(1,Method.Ncomponents); % Initializing the variable for the gradient of the components' wave functions in the z direction

%% Computation of the gradients
% FOR each component
for n = 1:Method.Ncomponents
        Grady_Phi{n} = ifft(FFTOperators3D.Gy.*fft(Phi{n},[],1),[],1); % First order derivating the component's wave function on the y direction via FFT and iFFT
        Gradx_Phi{n} = ifft(FFTOperators3D.Gx.*fft(Phi{n},[],2),[],2); % First order derivating the component's wave function on the x direction via FFT and iFFT
        Gradz_Phi{n} = ifft(FFTOperators3D.Gz.*fft(Phi{n},[],3),[],3); % First order derivating the component's wave function on the z direction via FFT and iFFT
end

%% Computation of the energy
% FOR each component
for n = 1:Method.Ncomponents
    %% Initialization 
    Transport_energy = zeros(FFTGeometry3D.Ny,FFTGeometry3D.Nx,FFTGeometry3D.Nz); % Initialization of the local transport energy
    Nonlinear_energy = zeros(FFTGeometry3D.Ny,FFTGeometry3D.Nx,FFTGeometry3D.Nz); % Initialization of the local nonlinear energy
    Potential_energy = zeros(FFTGeometry3D.Ny,FFTGeometry3D.Nx,FFTGeometry3D.Nz); % Initialization of the local potential energy
    Dispersion_energy = zeros(FFTGeometry3D.Ny,FFTGeometry3D.Nx,FFTGeometry3D.Nz); % Initialization of the local dispersion energy
    %% Computation of the transport, nonlinear and potential local energies
    % FOR each component where the gradient in the x direction is non null
    for m = FFTPhysics3D.Gradientx_function_Index{n}
        Transport_energy = Transport_energy + real(FFTPhysics3D.Gradientx_function{n,m}(FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z).*Gradx_Phi{m}.*conj(Phi{n})); % Computing the local potential energy for the gradient in the x direction for the component
    end
    % FOR each component where the gradient in the y direction is non null
    for m = FFTPhysics3D.Gradienty_function_Index{n}
        Transport_energy = Transport_energy + real(FFTPhysics3D.Gradienty_function{n,m}(FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z).*Grady_Phi{m}.*conj(Phi{n})); % Computing the local potential energy for the gradient in the y direction for the component
    end
    % FOR each component where the gradient in the z direction is non null
    for m = FFTPhysics3D.Gradientz_function_Index{n}
        Transport_energy = Transport_energy + real(FFTPhysics3D.Gradientz_function{n,m}(FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z).*Gradz_Phi{m}.*conj(Phi{n})); % Computing the local transport energy for the gradient in the z direction the component
    end
    % FOR each component where the potential is non null
    for m = FFTPhysics3D.Potential_function_Index{n}
        Potential_energy = Potential_energy + real(FFTPhysics3D.Potential_function{n,m}(FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z).*Phi{m}.*conj(Phi{n})); % Computing the local potential energy for the component
    end
    % FOR each component where the potential is non null
    for m = FFTPhysics3D.TimePotential_function_Index{n}
        Potential_energy = Potential_energy + real(FFTPhysics3D.TimePotential_function{n,m}(Method.Iterations*Method.Deltat,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z).*Phi{m}.*conj(Phi{n})); % Computing the local potential energy for the component
    end
    % FOR each component where the dispersion is non null
    for m = FFTPhysics3D.Dispersion_function_Index{n}
        Dispersion_energy = Dispersion_energy + real(ifftn(FFTPhysics3D.Dispersion{n,m}.*fftn(Phi{m})).*conj(Phi{n})); % Computing the wave function with the dispersion operator applied
    end
    % FOR each component where the dispersion is non null
    for m = FFTPhysics3D.TimeDispersion_function_Index{n}
        Dispersion_energy = Dispersion_energy + real(ifftn(FFTPhysics3D.TimeDispersion_function{n,m}(Method.Iterations*Method.Deltat,-1i*FFTOperators3D.Gx,-1i*FFTOperators3D.Gy,-1i*FFTOperators3D.Gz).*fftn(Phi{m})).*conj(Phi{n})); % Computing the wave function with the dispersion operator applied
    end
    % FOR each component where the stochastic potential is non null
    for m = FFTPhysics3D.StochasticPotential_function_Index{n}
        if (iscell(FFTPhysics3D.StochasticProcess_function))
            for m_noise = 1:length(FFTPhysics3D.StochasticProcess_function)
                FFTPhysics3D.StochasticProcess{m_noise} = (FFTPhysics3D.StochasticProcess_function{m_noise}((Method.Iterations + 1)*Method.Deltat,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z)-FFTPhysics3D.StochasticProcess_function{m_noise}(Method.Iterations*Method.Deltat,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z))/Method.Deltat;
            end
        else
            FFTPhysics3D.StochasticProcess = (FFTPhysics3D.StochasticProcess_function(Method.Iterations*Method.Deltat,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z) - FFTPhysics3D.StochasticProcess_function((Method.Iterations - 1)*Method.Deltat,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z))/Method.Deltat;
        end
        Potential_energy  = Potential_energy + real(FFTPhysics3D.StochasticPotential_function{n,m}(FFTPhysics3D.StochasticProcess,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z).*Phi{m}.*conj(Phi{n})); % Computing and storing the stochastic potential
    end
    % FOR each component where the nonlinearity is non null
    for m = FFTPhysics3D.Nonlinearity_function_Index{n}
        % IF the energy's nonlinearity is defined
        if (isempty(FFTPhysics3D.Nonlinearity_energy_function{n,m}) == 0)
            Nonlinear_energy = Nonlinear_energy + real(FFTPhysics3D.Beta*FFTPhysics3D.Nonlinearity_energy_function{n,m}(Phi,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z).*Phi{m}.*conj(Phi{n})); % Computing of the local nonlinear energy
        % ELSE if the energy's nonlinearity is not defined
        else
            Nonlinear_energy = Nonlinear_energy + real(FFTPhysics3D.Beta*FFTPhysics3D.Nonlinearity_function{n,m}(Phi,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z).*Phi{m}.*conj(Phi{n})); % Computing of the local nonlinear energy
        end
    end
    % FOR each component where the non-local nonlinearity is non null
    for m = FFTPhysics3D.FFTNonlinearity_function_Index{n}
        % IF the energy's non-local nonlinearity is defined
        if (isempty(FFTPhysics3D.FFTNonlinearity_energy_function{n,m}) == 0)
            Nonlinear_energy = Nonlinear_energy + real(FFTPhysics3D.FFTNonlinearity_energy_function{n,m}(Phi,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z,-1i*FFTOperators3D.Gx,-1i*FFTOperators3D.Gy,-1i*FFTOperators3D.Gz).*Phi{m}.*conj(Phi{n})); % Computing of the local non-local nonlinear energy
        % ELSE if the energy's non-local nonlinearity is not defined
        else
            Nonlinear_energy = Nonlinear_energy + real(FFTPhysics3D.FFTNonlinearity_function{n,m}(Phi,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z,-1i*FFTOperators3D.Gx,-1i*FFTOperators3D.Gy,-1i*FFTOperators3D.Gz).*Phi{m}.*conj(Phi{n})); % Computing of the local non-local nonlinear energy
        end
    end
    
    % FOR each component where the stochastic potential is non null
    for m = FFTPhysics3D.StochasticDispersion_function_Index{n}
        if (iscell(FFTPhysics3D.StochasticProcess_function))
            for m_noise = 1:length(FFTPhysics3D.StochasticProcess_function)
                FFTPhysics3D.StochasticProcess{m_noise} = (FFTPhysics3D.StochasticProcess_function{m_noise}((Method.Iterations + 1)*Method.Deltat,FFTGeometry1D.X,FFTGeometry1D.Y,FFTGeometry1D.Z)-FFTPhysics3D.StochasticProcess_function{m_noise}(Method.Iterations*Method.Deltat,FFTGeometry1D.X,FFTGeometry1D.Y,FFTGeometry1D.Z))/Method.Deltat;
            end
        else
            FFTPhysics3D.StochasticProcess = (FFTPhysics3D.StochasticProcess_function(Method.Iterations*Method.Deltat,FFTGeometry1D.X,FFTGeometry1D.Y,FFTGeometry1D.Z) - FFTPhysics3D.StochasticProcess_function((Method.Iterations - 1)*Method.Deltat,FFTGeometry1D.X,FFTGeometry1D.Y,FFTGeometry1D.Z))/Method.Deltat;
        end
        Dispersion_energy  = Dispersion_energy + real(ifftn(FFTPhysics3D.StochasticDispersion_function{n,m}(FFTPhysics3D.StochasticProcess,-1i*FFTOperators3D.Gx,-1i*FFTOperators3D.Gy,-1i*FFTOperators3D.Gz).*fftn(Phi{m})).*conj(Phi{n})); % Computing and storing the stochastic dispersion
    end
    %% Computation of the local energy
    Local_energy{n} = Potential_energy + Nonlinear_energy + Transport_energy + Dispersion_energy; % Computing local energy
    %% Integration over the space
    Energy{n} = FFTGeometry3D.dx*FFTGeometry3D.dy*FFTGeometry3D.dz*sum(sum(sum(Local_energy{n}))); % Integration of the local energy
end