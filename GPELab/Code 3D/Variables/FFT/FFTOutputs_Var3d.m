%% Computation of outputs and storage in the Outputs structure
%% INPUTS:
%%          FFTPhi: Wave functions in the 3D geometry for the FFT (cell array)
%%          Outputs: Different outputs computed during the computation of the ground states (structure) (see OutputsINI_Var3D.m)
%%          FFTGeometry3D: Structure containing variables concerning the geometry of the problem in 3D in the FFT context (structure) (see FFTGeometry3D_Var3d.m)
%%          FFTPhysics3D: Structure containing variables concerning the physics of the problem in 3D in the FFT context (structure) (see FFTPhysics3D_Var3d.m)
%%          FFTOperators3D: Structure containing the derivative FFT operators (structure) (see FFTOperators2D_Var3d.m)
%% OUTPUT:
%%          Outputs: Different outputs computed during the computation of the ground states (structure) (see OutputsINI_Var3D.m)
%% FUNCTIONS USED:
%%          Energy_GPE_Fourier3d: To compute the wave functions' energies using FFT operators (line 21)
%%          chemical_potential3d: To compute the wave function chemical potential (line 22)
%%          alpha_rms3d: To compute the root mean square in the x,y,z direction (line 26,27 and 28)
%%          Angular_momentum_Fourier3d: To compute the wave function angular momentum (line 31)

function [Outputs] = FFTOutputs_Var3d(FFTPhi, Outputs, Method, FFTGeometry3D, FFTPhysics3D, FFTOperators3D)
%% Updating iterations
Outputs.Iterations = Outputs.Iterations + 1;

%% Computing and storing outputs
Energy_tmp = Energy_GPE_Fourier3d(FFTPhi, Method, FFTGeometry3D, FFTPhysics3D, FFTOperators3D); % Computing the energy of the wave functions
Chemical_tmp = chemical_potential3d(FFTPhi, Energy_tmp, Method, FFTGeometry3D, FFTPhysics3D); % Computing the chemical potential of the wave functions
% IF there are user defined global functions
if (Outputs.User_compute_global)
    % FOR each user defined function
    for m = 1:Outputs.User_defined_number_global
        Outputs.User_defined_global{m}(Outputs.Iterations) = Outputs.User_defined_function_global{m}(FFTPhi,FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z,  -1i*FFTOperators3D.Gx, -1i*FFTOperators3D.Gy, -1i*FFTOperators3D.Gz); % Computing the user defined function applied to the wave function
    end
end
% FOR each component
for n = 1:Method.Ncomponents
    Outputs.phi_abs_0{n}(Outputs.Iterations)= abs(FFTPhi{n}(ceil((FFTGeometry3D.Ny+1)/2),ceil((FFTGeometry3D.Nx+1)/2),ceil((FFTGeometry3D.Nz+1)/2)))^2; % Computing of the square of the wave function at the origin
    Outputs.x_rms{n}(Outputs.Iterations) = alpha_rms3d(abs(FFTPhi{n}),FFTGeometry3D.X,FFTGeometry3D); % Computing of the rms in the x direction
    Outputs.y_rms{n}(Outputs.Iterations) = alpha_rms3d(abs(FFTPhi{n}),FFTGeometry3D.Y,FFTGeometry3D); % Computing of the rms in the y direction
    Outputs.z_rms{n}(Outputs.Iterations) = alpha_rms3d(abs(FFTPhi{n}),FFTGeometry3D.Z,FFTGeometry3D); % Computing of the rms in the z direction
    Outputs.Energy{n}(Outputs.Iterations) = Energy_tmp{n}; % Storing the energy of the wave function
    Outputs.Chemical_potential{n}(Outputs.Iterations) =  Chemical_tmp{n};  % Storing the chemical potential of the wave function
    Outputs.Angular_momentum{n}(Outputs.Iterations) = Angular_momentum_Fourier3d(FFTPhi{n}, FFTGeometry3D, FFTOperators3D); % Computing the angular momentum of the wave function
    % IF one has chosen to save the functions during the simulation
    if (Outputs.Save_solution)
        Phi = zeros(FFTGeometry3D.Ny+1,FFTGeometry3D.Nx+1,FFTGeometry3D.Nz+1); % Initializing the variable for the storing of the function
        Phi(1:FFTGeometry3D.Ny,1:FFTGeometry3D.Nx,1:FFTGeometry3D.Nz) = FFTPhi{n}; % Storing of the ground states solutions
        Phi(FFTGeometry3D.Ny+1,1:FFTGeometry3D.Nx,1:FFTGeometry3D.Nz) = FFTPhi{n}(1,:,:);
        Phi(1:FFTGeometry3D.Ny,FFTGeometry3D.Nx+1,1:FFTGeometry3D.Nz) = FFTPhi{n}(:,1,:);
        Phi(1:FFTGeometry3D.Ny,1:FFTGeometry3D.Nx,FFTGeometry3D.Nz+1) = FFTPhi{n}(:,:,1);
        Phi(FFTGeometry3D.Ny+1,FFTGeometry3D.Nx+1,FFTGeometry3D.Nz+1) = FFTPhi{n}(1,1,1);
        Outputs.Solution{Outputs.Iterations}{n} = Phi; % Saving the function in the outputs
    end
    % IF there are user defined local functions
    if (Outputs.User_compute_local)
        % FOR each user defined function
        for m = 1:Outputs.User_defined_number_local
            Outputs.User_defined_local{n,m}(Outputs.Iterations) = Outputs.User_defined_function_local{m}(FFTPhi{n},FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z,  -1i*FFTOperators3D.Gx, -1i*FFTOperators3D.Gy, -1i*FFTOperators3D.Gz); % Computing the user defined function applied to the wave function
        end
    end
end