%% Computation of outputs and storage in the Outputs structure
%% INPUTS:
%%          FDPhi: Wave functions in the 1D geometry for the FD (cell array)
%%          Outputs: Different outputs computed during the computation of the ground states (structure) (see OutputsINI_Var1D.m)
%%          FDGeometry1D: Structure containing variables concerning the geometry of the problem in 1D in the FD context (structure) (see FDGeometry1D_Var1d.m)
%%          FDPhysics1D: Structure containing variables concerning the physics of the problem in 1D in the FD context (structure) (see FDPhysics1D_Var1d.m)
%%          FDOperators1D: Structure containing the derivative FFT operators (structure) (see FDOperators1D_Var1d.m)
%% OUTPUT:
%%          Outputs: Different outputs computed during the computation of the ground states (structure) (see OutputsINI_Var1D.m)
%% FUNCTIONS USED:
%%          Energy_GPE_Diff1d: To compute the wave functions' energies using FFT operators (line 20)
%%          chemical_potential1d: To compute the wave function chemical potential (line 21)
%%          alpha_rms1d: To compute the root mean square in the x direction (line 26)

function [Outputs] = FDOutputs_Var1d(FDPhi, Outputs, Method, FDGeometry1D, FDPhysics1D, FDOperators1D)
%% Updating iterations
Outputs.Iterations = Outputs.Iterations +1;

%% Computing and storing outputs
Energy_tmp = Energy_GPE_Diff1d(FDPhi, Method, FDGeometry1D, FDPhysics1D, FDOperators1D); % Computing the energy of the wave functions
Chemical_tmp = Chemical_potential1d(FDPhi, Energy_tmp, Method, FDGeometry1D, FDPhysics1D); % Computing the chemical potential of the wave functions
% IF there are user defined functions
if (Outputs.User_compute_global)
    % FOR each user defined function
    for m = 1:Outputs.User_defined_number_global
        Outputs.User_defined_global{m}(Outputs.Iterations) = Outputs.User_defined_function_global{m}(FDPhi,FDGeometry1D.X); % Computing the user defined function applied to the wave function
    end
end
% FOR each component
for n = 1:Method.Ncomponents
    Outputs.phi_abs_0{n}(Outputs.Iterations)= abs(FDPhi{n}(ceil((FDGeometry1D.Nx+1)/2)))^2; % Computing of the square of the wave function at the origin
    Outputs.x_rms{n}(Outputs.Iterations) = alpha_rms1d(abs(FDPhi{n}),FDGeometry1D.X,FDGeometry1D); % Computing of the rms in the x direction
    Outputs.Energy{n}(Outputs.Iterations) = Energy_tmp{n}; % Storing the energy of the wave function
    Outputs.Chemical_potential{n}(Outputs.Iterations) =  Chemical_tmp{n}; % Storing the chemical potential of the wave function
    % IF one has chosen to save the functions during the simulation
    if (Outputs.Save_solution)
        Phi = zeros(1,FDGeometry2D.Nx+2); % Initializing the variable for the storing of the function
        Phi(2:FDGeometry2D.Nx) = FDPhi{n}; % Storing of the function
        Outputs.Solution{Outputs.Iterations}{n} = Phi; % Saving the function in the outputs
    end
    % IF there are user defined functions
    if (Outputs.User_compute_local)
        % FOR each user defined function
        for m = 1:Outputs.User_defined_number_local
            Outputs.User_defined_local{n,m}(Outputs.Iterations) = Outputs.User_defined_function_local{m}(FDPhi{n},FDGeometry1D.X); % Computing the user defined function applied to the wave function
        end
    end
end