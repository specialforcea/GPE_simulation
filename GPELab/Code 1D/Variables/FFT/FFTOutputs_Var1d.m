%% Computation of outputs and storage in the Outputs structure
%% INPUTS:
%%          FFTPhi: Wave functions in the 3D geometry for the FFT (cell array)
%%          Outputs: Different outputs computed during the computation of the ground states (structure) (see OutputsINI_Var1D.m)
%%          FFTGeometry1D: Structure containing variables concerning the geometry of the problem in 1D in the FFT context (structure) (see FFTGeometry1D_Var1d.m)
%%          FFTPhysics1D: Structure containing variables concerning the physics of the problem in 1D in the FFT context (structure) (see FFTPhysics1D_Var1d.m)
%%          FFTOperators1D: Structure containing the derivative FFT operators (structure) (see FFTOperators1D_Var1d.m)
%% OUTPUT:
%%          Outputs: Different outputs computed during the computation of the ground states (structure) (see OutputsINI_Var1D.m)
%% FUNCTIONS USED:
%%          Energy_GPE_Fourier1d: To compute the wave functions' energies using FFT operators (line 20)
%%          chemical_potential1d: To compute the wave function chemical potential (line 21)
%%          alpha_rms1d: To compute the root mean square in the x direction (line 25)

function [Outputs] = FFTOutputs_Var1d(FFTPhi, Outputs, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D)
%% Updating iterations
Outputs.Iterations = Outputs.Iterations + 1;

%% Computing and storing outputs
Energy_tmp = Energy_GPE_Fourier1d(FFTPhi, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D); % Computing the energy of the wave functions
Chemical_tmp = Chemical_potential1d(FFTPhi, Energy_tmp, Method, FFTGeometry1D, FFTPhysics1D); % Computing the chemical potential of the wave functions
% IF there are user defined functions
if (Outputs.User_compute_global)
    % FOR each user defined function
    for m = 1:Outputs.User_defined_number_global
        Outputs.User_defined_global{m}(Outputs.Iterations) = Outputs.User_defined_function_global{m}(FFTPhi,FFTGeometry1D.X, -1i*FFTOperators1D.Gx); % Computing the user defined function applied to the wave function
    end
end
% FOR each component
for n = 1:Method.Ncomponents
    Outputs.phi_abs_0{n}(Outputs.Iterations)= abs(FFTPhi{n}(ceil((FFTGeometry1D.Nx+1)/2)))^2; % Computing of the square of the wave function at the origin
    Outputs.x_rms{n}(Outputs.Iterations) = alpha_rms1d(abs(FFTPhi{n}),FFTGeometry1D.X,FFTGeometry1D); % Computing of the rms in the x direction
    Outputs.Energy{n}(Outputs.Iterations) = Energy_tmp{n}; % Storing the energy of the wave function
    Outputs.Chemical_potential{n}(Outputs.Iterations) =  Chemical_tmp{n};  % Storing the chemical potential of the wave function
    % IF one has chosen to save the functions during the simulation
    if (Outputs.Save_solution)
        Phi = zeros(1,FFTGeometry1D.Nx+1); % Initializing the variable for the storing of the function
        Phi(1:FFTGeometry1D.Nx) = FFTPhi{n}; % Storing of the function
        Outputs.Solution{Outputs.Iterations}{n} = Phi; % Saving the function in the outputs
    end
    % IF there are user defined functions
    if (Outputs.User_compute_local)
        % FOR each user defined function
        for m = 1:Outputs.User_defined_number_local
            Outputs.User_defined_local{n,m}(Outputs.Iterations) = Outputs.User_defined_function_local{m}(FFTPhi{n},FFTGeometry1D.X, -1i*FFTOperators1D.Gx); % Computing the user defined function applied to the wave function
        end
    end
end