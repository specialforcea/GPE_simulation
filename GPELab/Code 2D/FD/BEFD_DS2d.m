%% Backward Euler Finite Differences for the Dynamic System (BEFD-DS) : Semi-implicit FD scheme for the computation of the dynamic Gross Pitaevskii equation
%% INPUTS:
%%          Phi_0: Initial wave functions (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          Geometry2D: Structure containing variables concerning the geometry of the problem in 2D (structure) (see Geometry2D_Var2d.m)
%%          Physics2D: Structure containing variables concerning the physics of the problem in 2D (structure) (see Physics2D_Var2d.m)
%%          Outputs: Different outputs computed during the computation of the problem (structure) (see OutputsINI_Var2D.m and FFTOutputs_Var2d.m)
%%          Print: Structure containing variables concerning the printing and drawing of informations during the computations (structure) (see Print_Var2d.m)
%%          Figure: Structure containing variables concerning the figures (structure)
%% OUTPUTS:
%%          Phi_final: Dynamic's wave functions computated with the BEFD method (cell array)
%%          Outputs: Different outputs computed during the computation of the problem (structure) (see OutputsINI_Var2D.m and FFTOutputs_Var2d.m)
%% FUNCTIONS USED:
%%          FDGeometry2D_Var2d: To compute the 2d geometry used for the FD (line 30)
%%          FDPhysics2D_Var2d: To compute the 2d physics used for the FD (line 31)
%%          FDOperators2D_Var2d: To compute the 2d operators used for the FD (line 39)
%%          Local_Krylov_BEFD_solution2d: To compute the CNFG on a single time step using the iterative BEFD method (line 57)
%%          Local_Direct_BEFD_solution2d: To compute the CNFG on a single time step using the direct BEFD method (line 63)
%%          FDOutputs_Var2d: To compute output variables (line 74,78 and 91)
%%          Print_Info2d: To print informations concerning the wave functions (line 81 and 92)
%%          Draw_solution2d: To draw the ground states (line 85)

function [Phi_final, Outputs] = BEFD_DS2d(Phi_0, Method, Geometry2D, Physics2D, Outputs, Print, Figure)

%% CPUtime initilization
Method.Cputime = 0; % Initialization of the CPUtime variable
Method.Cputime_temp = cputime; % Initialization of the relative CPUtime variable

%% Geometry, physics and ground states initialization for FD
FDGeometry2D = FDGeometry2D_Var2d(Geometry2D); % Changing the geometry for the FD
FDPhysics2D = FDPhysics2D_Var2d(Method, Physics2D, FDGeometry2D); % Changing the physics for the FD
% Changing the solution for the FD
% FOR each component
for n = 1:Method.Ncomponents
    FDPhi{n} = Phi_0{n}(2:FDGeometry2D.Ny+1,2:FDGeometry2D.Nx+1); % Removing boundaries for the FD
end

%% Initialization of the derivative FD operators
FDOperators2D = FDOperators2D_Var2d(Method, FDGeometry2D, FDPhysics2D); % Computing the derivative FD operators

%% CPUtime initilization bis 
Method.Cputime_temp = cputime - Method.Cputime_temp; % Storing the CPUtime relatively to the begining of the program
Method.Cputime  = Method.Cputime + Method.Cputime_temp;  % Updating the CPUtime variable

%% Computation of dynamic of the GPE via BEFD
% Stopping criterions: stopping time
while (Method.Iterations*Method.Deltat<Method.Stop_time)
    %% Updating variables
    FDPhi_tmp = FDPhi; % Storing a temporary variable of the ground states to compute local evolution
    Method.Cputime_temp = cputime; % Reinitialization of the relative CPUtime variable
    Method.Iterations = Method.Iterations + 1; % Iteration count

    %% BEFD Type: iterative method or direct inversion of the system
    % IF one has chosen to use the iterative method
    if (Method.Solver_FD == 0)
    %% Computing ground state using the BEFD method with iterative method
    [FDPhi, flag, relres, iter, resvec] = Local_Krylov_BEFD_solution2d(FDPhi, Method, FDGeometry2D, FDPhysics2D, FDOperators2D); % Computation of the ground state using the iterative BEFD method on a single step of time
    Outputs.iteration_vec(Method.Iterations) = size(resvec,1); % Storing the total number of BICGSTAB iterations at the iteration-th iteration of the Gradient Flow

    % ELSEIF one has chosen to use the direct inversion
    elseif (Method.Solver_FD == 1)
    %% Computing ground state using the BEFD method with direct inversion
    FDPhi = Local_Direct_BEFD_solution2d(FDPhi, Method, FDGeometry2D, FDPhysics2D, FDOperators2D); % Computation of the ground state using the direct BEFD method on a single step of time
    end
    
    %% Updating CPUtime
    Method.Cputime_temp = cputime - Method.Cputime_temp; % Storing the CPUtime relatively to the begining of the "while" loop
    Method.Cputime = Method.Cputime + Method.Cputime_temp; % Updating the CPUtime variable
    
    %% Computation of outputs
    % IF one wants to either compute outputs or to print informations during 
    % the computation
    if (Method.Output) && (mod(Method.Iterations,Outputs.Evo_outputs) == 0)
        Outputs = FDOutputs_Var2d(FDPhi, Outputs, Method, FDGeometry2D, FDPhysics2D, FDOperators2D); % Computing the output variables
    end
    %% Printing informations and drawing ground states
    % IF the number of iterations has been reached
    if (mod(Method.Iterations,Print.Evo) == 0)
        % IF one has chosen to print informations
        if (Method.Output) && (Print.Print == 1)
            Print_Info2d(Outputs, Method) % Printing informations
        elseif (~Method.Output) && (Print.Print == 1)
            Outputstmp = FDOutputs_Var2d(FDPhi, Outputs, Method, FDGeometry2D, FDPhysics2D, FDOperators2D); % Computing the output variables
            Print_Info2d(Outputstmp, Method) % Printing informations
        end
        % IF one has chosen to draw the ground states
        if (Print.Draw == 1)
            Draw_solution2d(FDPhi,Method,FDGeometry2D,Figure) % Drawing the wave functions' modulus square and angle
        end;
    end
end

%% Printing the informations of the final ground states
Outputs = FDOutputs_Var2d(FDPhi, Outputs, Method, FDGeometry2D, FDPhysics2D, FDOperators2D); % Computing the output variables
Print_Info2d(Outputs, Method) % Printing informations

%% Finalization of the ground states
% FOR each component
for n = 1:Method.Ncomponents
Phi_final{n} = zeros(Geometry2D.Ny,Geometry2D.Nx); % Initializing the variable for the storing of the final ground states
Phi_final{n}(2:FDGeometry2D.Ny+1,2:FDGeometry2D.Nx+1) = FDPhi{n}; % Storing of the ground states solutions
end