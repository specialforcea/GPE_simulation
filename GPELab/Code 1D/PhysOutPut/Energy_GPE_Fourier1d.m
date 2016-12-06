%% Computation of the energy of the GPE using the FFT operators
%% INPUTS:
%%          Phi: Wave functions (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var1d.m)
%%          FFTGeometry1D: Structure containing variables concerning the geometry of the problem in 1D in the FFT context (structure) (see FFTGeometry1D_Var1d.m)
%%          FFTPhysics1D: Structure containing variables concerning the physics of the problem in 1D in the FFT context (structure) (see FFTPhysics1D_Var1d.m)
%%          FFTOperators1D: Structure containing the derivative FFT operators (structure) (see FFTOperators1D_Var1d.m)
%% OUTPUT:
%%          Energy: Energy of the wave functions (cell array)
%% FUNCTIONS USED:
%%          L2_norm1d: To integrate the local energy (line 62)

function [Energy] = Energy_GPE_Fourier1d(Phi, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D)
%% Initialization
Gradx_Phi = cell(1,Method.Ncomponents); % Initializing the variable for the gradient of the components' wave functions in the x direction

%% Computation of the gradients
% FOR each component
for n = 1:Method.Ncomponents
        Gradx_Phi{n} = ifft(FFTOperators1D.Gx.*fft(Phi{n})); % First order derivating the component's wave function on the x direction via FFT and iFFT
end

%% Computation of the energy
% FOR each component
for n = 1:Method.Ncomponents
    %% Initialization 
    Transport_energy = zeros(FFTGeometry1D.Nx,1); % Initialization of the local transport energy
    Nonlinear_energy = zeros(FFTGeometry1D.Nx,1); % Initialization of the local nonlinear energy
    Potential_energy = zeros(FFTGeometry1D.Nx,1); % Initialization of the local potential energy
    Dispersion_energy = zeros(FFTGeometry1D.Nx,1); % Initialization of the local dispersion energy
    %% Computation of the transport, nonlinear and potential local energies
    % FOR each component where the gradient in the x direction is non null
    for m = FFTPhysics1D.Gradientx_function_Index{n}
        Transport_energy = Transport_energy + real(FFTPhysics1D.Gradientx_function{n,m}(FFTGeometry1D.X).*Gradx_Phi{m}.*conj(Phi{n})); % Computing the local potential energy for the gradient in the x direction for the component
    end
    % FOR each component where the potential is non null
    for m = FFTPhysics1D.Potential_function_Index{n}
        Potential_energy = Potential_energy + real(FFTPhysics1D.Potential_function{n,m}(FFTGeometry1D.X).*Phi{m}.*conj(Phi{n})); % Computing the local potential energy for the component
    end
    % FOR each component where the time-dependent potential is non null
    for m = FFTPhysics1D.TimePotential_function_Index{n}
        Potential_energy = Potential_energy + real(FFTPhysics1D.TimePotential_function{n,m}(Method.Iterations*Method.Deltat,FFTGeometry1D.X).*Phi{m}.*conj(Phi{n})); % Computing the local time-dependent potential energy for the component
    end
    % FOR each component where the stochastic potential is non null
    for m = FFTPhysics1D.StochasticPotential_function_Index{n}
        if (iscell(FFTPhysics1D.StochasticProcess_function))
            for m_noise = 1:length(FFTPhysics1D.StochasticProcess_function)
                FFTPhysics1D.StochasticProcess{m_noise} = (FFTPhysics1D.StochasticProcess_function{m_noise}(Method.Iterations*Method.Deltat,FFTGeometry1D.X)-FFTPhysics1D.StochasticProcess_function{m_noise}((Method.Iterations - 1)*Method.Deltat,FFTGeometry1D.X))/Method.Deltat;
            end
        else
            FFTPhysics1D.StochasticProcess = (FFTPhysics1D.StochasticProcess_function(Method.Iterations*Method.Deltat,FFTGeometry1D.X) - FFTPhysics1D.StochasticProcess_function((Method.Iterations - 1)*Method.Deltat,FFTGeometry1D.X))/Method.Deltat;
        end
        Potential_energy  = Potential_energy + real(FFTPhysics1D.StochasticPotential_function{n,m}(FFTPhysics1D.StochasticProcess,FFTGeometry1D.X).*Phi{m}.*conj(Phi{n})); % Computing and storing the stochastic potential
    end
    % FOR each component where the nonlinearity is non null
    for m = FFTPhysics1D.Nonlinearity_function_Index{n}
        % IF the energy's nonlinearity is defined
        if (isempty(FFTPhysics1D.Nonlinearity_energy_function{n,m}) == 0)
            Nonlinear_energy = Nonlinear_energy + real(FFTPhysics1D.Beta*FFTPhysics1D.Nonlinearity_energy_function{n,m}(Phi,FFTGeometry1D.X).*Phi{m}.*conj(Phi{n})); % Computing of the local nonlinear energy
        % ELSE if the energy's nonlinearity is not defined
        else
            Nonlinear_energy = Nonlinear_energy + real(FFTPhysics1D.Beta*FFTPhysics1D.Nonlinearity_function{n,m}(Phi,FFTGeometry1D.X).*Phi{m}.*conj(Phi{n})); % Computing of the local nonlinear energy
        end
    end
    % FOR each component where the non-local nonlinearity is non null
    for m = FFTPhysics1D.FFTNonlinearity_function_Index{n}
        % IF the energy's non-local nonlinearity is defined
        if (isempty(FFTPhysics1D.FFTNonlinearity_energy_function{n,m}) == 0)
            Nonlinear_energy = Nonlinear_energy + real(FFTPhysics1D.FFTNonlinearity_energy_function{n,m}(Phi,FFTGeometry1D.X,-1i*FFTOperators1D.Gx).*Phi{m}.*conj(Phi{n})); % Computing of the local non-local nonlinear energy
        % ELSE if the energy's non-local nonlinearity is not defined
        else
            Nonlinear_energy = Nonlinear_energy + real(FFTPhysics1D.FFTNonlinearity_function{n,m}(Phi,FFTGeometry1D.X,-1i*FFTOperators1D.Gx).*Phi{m}.*conj(Phi{n})); % Computing of the local non-local nonlinear energy
        end
    end
    % FOR each component where the dispersion is non null
    for m = FFTPhysics1D.Dispersion_function_Index{n}
        Dispersion_energy  = Dispersion_energy + real(ifft(FFTPhysics1D.Dispersion_function{n,m}(-1i*FFTOperators1D.Gx).*fft(Phi{m})).*conj(Phi{n})); % Computing and storing the dispersion
    end
    % FOR each component where the dispersion is non null
    for m = FFTPhysics1D.TimeDispersion_function_Index{n}
        Dispersion_energy  = Dispersion_energy + real(ifft(FFTPhysics1D.TimeDispersion_function{n,m}(Method.Iterations*Method.Deltat,-1i*FFTOperators1D.Gx).*fft(Phi{m})).*conj(Phi{n})); % Computing and storing the dispersion
    end
    % FOR each component where the stochastic dispersion is non null
    for m = FFTPhysics1D.StochasticDispersion_function_Index{n}
        if (iscell(FFTPhysics1D.StochasticProcess_function))
            for m_noise = 1:length(FFTPhysics1D.StochasticProcess_function)
                FFTPhysics1D.StochasticProcess{m_noise} = (FFTPhysics1D.StochasticProcess_function{m_noise}((Method.Iterations + 1)*Method.Deltat,FFTGeometry1D.X)-FFTPhysics1D.StochasticProcess_function{m_noise}(Method.Iterations*Method.Deltat,FFTGeometry1D.X))/Method.Deltat;
            end
        else
            FFTPhysics1D.StochasticProcess = (FFTPhysics1D.StochasticProcess_function(Method.Iterations*Method.Deltat,FFTGeometry1D.X) - FFTPhysics1D.StochasticProcess_function((Method.Iterations - 1)*Method.Deltat,FFTGeometry1D.X))/Method.Deltat;
        end
        Dispersion_energy  = Dispersion_energy + real(ifft(FFTPhysics1D.StochasticDispersion_function{n,m}(FFTPhysics1D.StochasticProcess,-1i*FFTOperators1D.Gx).*fft(Phi{m})).*conj(Phi{n})); % Computing and storing the stochastic potential
    end
    %% Computation of the local energy
    Local_energy{n} = Potential_energy + Nonlinear_energy + Transport_energy + Dispersion_energy ; % Computing local energy
    %% Integration over the space
    Energy{n} = FFTGeometry1D.dx*sum(Local_energy{n}); % Integration of the local energy
end