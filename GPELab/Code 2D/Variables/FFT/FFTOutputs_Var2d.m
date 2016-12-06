%% Computation of outputs and storage in the Outputs structure
%% INPUTS:
%%          FFTPhi: Wave functions in the 3D geometry for the FFT (cell array)
%%          Outputs: Different outputs computed during the computation of the ground states (structure) (see OutputsINI_Var2D.m)
%%          FFTGeometry2D: Structure containing variables concerning the geometry of the problem in 2D in the FFT context (structure) (see FFTGeometry2D_Var2d.m)
%%          FFTPhysics2D: Structure containing variables concerning the physics of the problem in 2D in the FFT context (structure) (see FFTPhysics2D_Var2d.m)
%%          FFTOperators2D: Structure containing the derivative FFT operators (structure) (see FFTOperators2D_Var2d.m)
%% OUTPUT:
%%          Outputs: Different outputs computed during the computation of the ground states (structure) (see OutputsINI_Var2D.m)
%% FUNCTIONS USED:
%%          Energy_GPE_Fourier2d: To compute the wave functions' energies using FFT operators (line 21)
%%          chemical_potential2d: To compute the wave function chemical potential (line 22)
%%          alpha_rms2d: To compute the root mean square in the x,y direction (line 26 and 27)
%%          Angular_momentum_Fourier2d: To compute the wave function angular momentum (line 30)

function [Outputs] = FFTOutputs_Var2d(FFTPhi, Outputs, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D)
%% Updating iterations
Outputs.Iterations = Outputs.Iterations + 1;

%% Computing and storing outputs
Energy_tmp = Energy_GPE_Fourier2d(FFTPhi, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D); % Computing the energy of the wave functions
Chemical_tmp = chemical_potential2d(FFTPhi, Energy_tmp, Method, FFTGeometry2D, FFTPhysics2D); % Computing the chemical potential of the wave functions
% IF there are user defined global functions
if (Outputs.User_compute_global)
    % FOR each user defined function
    for m = 1:Outputs.User_defined_number_global
        Outputs.User_defined_global{m}(Outputs.Iterations) = Outputs.User_defined_function_global{m}(FFTPhi,FFTGeometry2D.X,FFTGeometry2D.Y, -1i*FFTOperators2D.Gx, -1i*FFTOperators2D.Gy); % Computing the user defined function applied to the wave function
    end
end
% FOR each component
for n = 1:Method.Ncomponents
    Outputs.phi_abs_0{n}(Outputs.Iterations)= abs(FFTPhi{n}(ceil((FFTGeometry2D.Ny+1)/2),ceil((FFTGeometry2D.Nx+1)/2)))^2; % Computing of the square of the wave function at the origin
    Outputs.x_rms{n}(Outputs.Iterations) = alpha_rms2d(abs(FFTPhi{n}),FFTGeometry2D.X,FFTGeometry2D); % Computing of the rms in the x direction
    Outputs.y_rms{n}(Outputs.Iterations) = alpha_rms2d(abs(FFTPhi{n}),FFTGeometry2D.Y,FFTGeometry2D); % Computing of the rms in the y direction
    Outputs.Energy{n}(Outputs.Iterations) = Energy_tmp{n}; % Storing the energy of the wave function
    Outputs.Chemical_potential{n}(Outputs.Iterations) =  Chemical_tmp{n};  % Storing the chemical potential of the wave function
    Outputs.Angular_momentum{n}(Outputs.Iterations) = Angular_momentum_Fourier2d(FFTPhi{n}, FFTGeometry2D, FFTOperators2D); % Computing the angular momentum of the wave function
    % IF one has chosen to save the functions during the simulation
    if (Outputs.Save_solution)
        Phi = zeros(FFTGeometry2D.Ny+1,FFTGeometry2D.Nx+1); % Initializing the variable for the storing of the function
        Phi(1:FFTGeometry2D.Ny,1:FFTGeometry2D.Nx) = FFTPhi{n}; % Storing of the function
        Outputs.Solution{Outputs.Iterations}{n} = Phi; % Saving the function in the outputs
    end
    % IF there are user defined local functions
    if (Outputs.User_compute_local)
        % FOR each user defined function
        for m = 1:Outputs.User_defined_number_local
            Outputs.User_defined_local{n,m}(Outputs.Iterations) = Outputs.User_defined_function_local{m}(FFTPhi{n},FFTGeometry2D.X,FFTGeometry2D.Y, -1i*FFTOperators2D.Gx, -1i*FFTOperators2D.Gy); % Computing the user defined function applied to the wave function
        end
    end
end
    