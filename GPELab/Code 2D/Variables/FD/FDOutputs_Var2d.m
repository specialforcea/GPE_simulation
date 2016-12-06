%% Computation of outputs and storage in the Outputs structure
%% INPUTS:
%%          FDPhi: Wave functions in the 2D geometry for the FD (cell array)
%%          Outputs: Different outputs computed during the computation of the ground states (structure) (see OutputsINI_Var2D.m)
%%          FDGeometry2D: Structure containing variables concerning the geometry of the problem in 2D in the FD context (structure) (see FDGeometry2D_Var2d.m)
%%          FDPhysics2D: Structure containing variables concerning the physics of the problem in 2D in the FD context (structure) (see FDPhysics2D_Var2d.m)
%%          FDOperators2D: Structure containing the derivative FFT operators (structure) (see FDOperators2D_Var2d.m)
%% OUTPUT:
%%          Outputs: Different outputs computed during the computation of the ground states (structure) (see OutputsINI_Var2D.m)
%% FUNCTIONS USED:
%%          Energy_GPE_Diff2d: To compute the wave functions' energies using FFT operators (line 20)
%%          chemical_potential2d: To compute the wave function chemical potential (line 21)
%%          Angular_momentum_Diff2d: To compute the wave function angular momentum (line 22)
%%          alpha_rms2d: To compute the root mean square in the x,y direction (line 27 and 28)

function [Outputs] = FDOutputs_Var2d(FDPhi, Outputs, Method, FDGeometry2D, FDPhysics2D, FDOperators2D)
%% Updating iterations
Outputs.Iterations = Outputs.Iterations + 1;
%% Computing and storing outputs
Energy_tmp = Energy_GPE_Diff2d(FDPhi, Method, FDGeometry2D, FDPhysics2D, FDOperators2D); % Computing the energy of the wave functions
Chemical_tmp = chemical_potential2d(FDPhi, Energy_tmp, Method, FDGeometry2D, FDPhysics2D); % Computing the chemical potential of the wave functions
Angular_tmp = Angular_momentum_Diff2d(FDPhi, Method, FDGeometry2D, FDOperators2D); % Computing the angular momentum of the wave function
% IF there are user defined local functions
if (Outputs.User_compute_global)
    % FOR each user defined function
    for m = 1:Outputs.User_defined_number_global
        Outputs.User_defined_global{m}(Outputs.Iterations) = Outputs.User_defined_function_global{m}(FDPhi,FDGeometry2D.X,FDGeometry2D.Y); % Computing the user defined function applied to the wave function
    end
end
% FOR each component
for n = 1:Method.Ncomponents
    Outputs.phi_abs_0{n}(Outputs.Iterations)= abs(FDPhi{n}(ceil((FDGeometry2D.Ny+1)/2),ceil((FDGeometry2D.Nx+1)/2)))^2; % Computing of the square of the wave function at the origin
    Outputs.x_rms{n}(Outputs.Iterations) = alpha_rms2d(abs(FDPhi{n}),FDGeometry2D.X,FDGeometry2D); % Computing of the rms in the x direction
    Outputs.y_rms{n}(Outputs.Iterations) = alpha_rms2d(abs(FDPhi{n}),FDGeometry2D.Y,FDGeometry2D); % Computing of the rms in the y direction
    Outputs.Energy{n}(Outputs.Iterations) = Energy_tmp{n}; % Storing the energy of the wave function
    Outputs.Chemical_potential{n}(Outputs.Iterations) =  Chemical_tmp{n}; % Storing the chemical potential of the wave function
    Outputs.Angular_momentum{n}(Outputs.Iterations) = Angular_tmp{n}; % Storing the angular momentum of the wave function
    % IF one has chosen to save the functions during the simulation
    if (Outputs.Save_solution)
        Phi = zeros(FDGeometry2D.Ny+2,FDGeometry2D.Nx+2); % Initializing the variable for the storing of the function
        Phi(2:FDGeometry2D.Ny+1,2:FDGeometry2D.Nx+1) = FDPhi{n}; % Storing of the function
        Outputs.Solution{Outputs.Iterations}{n} = Phi; % Saving the function in the outputs
    end
    % IF there are user defined local functions
    if (Outputs.User_compute_local)
        % FOR each user defined function
        for m = 1:Outputs.User_defined_number_global
            Outputs.User_defined_local{n,m}(Outputs.Iterations) = Outputs.User_defined_function_local{m}(FDPhi{n},FDGeometry2D.X,FDGeometry2D.Y); % Computing the user defined function applied to the wave function
        end
    end
end