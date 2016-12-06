%% Cranl Nicholson Finite Differences for the Continuous Normalized Gradient Flow (CNFD-CNFG) : Semi-implicit FD scheme for the computation of grounds states of the Gross Pitaevskii equation
%% INPUTS:
%%          Phi_0: Initial wave functions (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var1d.m)
%%          Geometry1D: Structure containing variables concerning the geometry of the problem in 1D (structure) (see Geometry1D_Var1d.m)
%%          Physics1D: Structure containing variables concerning the physics of the problem in 1D (structure) (see Physics1D_Var1d.m)
%%          Outputs: Different outputs computed during the computation of the ground states (structure) (see OutputsINI_Var1D.m and FFTOutputs_Var1d.m)
%%          Print: Structure containing variables concerning the printing and drawing of informations during the computations (structure) (see Print_Var1d.m)
%%          Figure: Structure containing variables concerning the figures (structure)
%% OUTPUTS:
%%          Phi_final: Ground states' wave functions computated with the CNFD method (cell array)
%%          Outputs: Different outputs computed during the computation of the ground states (structure) (see OutputsINI_Var1D.m and FFTOutputs_Var1d.m)
%% FUNCTIONS USED:
%%          FDGeometry1D_Var1d: To compute the 1d geometry used for the FD (line 31)
%%          FDPhysics1D_Var1d: To compute the 1d physics used for the FD (line 32)
%%          FDOperators1D_Var1d: To compute the 1d operators used for the FD (line 40)
%%          Local_Krylov_CNFD_solution1d: To compute the CNFG on a single time step using the iterative CNFD method (line 63)
%%          Local_Direct_CNFD_solution1d: To compute the CNFG on a single time step using the direct CNFD method (line 69)
%%          L2_norm1d: To compute the normalized wave functions (line 75)
%%          FDOutputs_Var1d: To compute output variables (line 88,92 and 105)
%%          Print_Info1d: To print informations concerning the wave functions (line 95 and 106)
%%          Draw_solution1d: To draw the ground states (line 99)

function [Phi_final, Outputs] = CNFD_CNGF1d(Phi_0, Method, Geometry1D, Physics1D, Outputs, Print, Figure)

%% CPUtime initilization
Method.Cputime = 0; % Initialization of the CPUtime variable
Method.Cputime_temp = cputime; % Initialization of the relative CPUtime variable

%% Geometry, physics and ground states initialization for FD
FDGeometry1D = FDGeometry1D_Var1d(Geometry1D); % Changing the geometry for the FD
FDPhysics1D = FDPhysics1D_Var1d(Method, Physics1D, FDGeometry1D); % Changing the physics for the FD
% Changing the solution for the FD
% FOR each component
for n = 1:Method.Ncomponents
    FDPhi{n} = Phi_0{n}(2:FDGeometry1D.Nx+1); % Removing boundaries for the FD
end

%% Initialization of the derivative FD operators
FDOperators1D = FDOperators1D_Var1d(Method, FDGeometry1D, FDPhysics1D); % Computing the derivative FD operators

%% CPUtime initilization bis 
Method.Cputime_temp = cputime - Method.Cputime_temp; % Storing the CPUtime relatively to the begining of the program
Method.Cputime  = Method.Cputime + Method.Cputime_temp;  % Updating the CPUtime variable

%% Computation of the ground state via CNFD
% Stopping criterions:- the global evolution of the ground states (maximum
%                     of each local evolution) must be lower than the 
%                     stopping criterion Method.Stop multiplied by the
%                     time step
%                     - the number of iterations must not exceed the
%                     maximum number of iterations  
while (Method.EvolutionCriterion>Method.Stop_crit{2}*Method.Deltat)&&(Method.Iterations<Method.Max_iter)
    %% Updating variables
    FDPhi_tmp = FDPhi; % Storing a temporary variable of the ground states to compute local evolution
    Method.Cputime_temp = cputime; % Reinitialization of the relative CPUtime variable
    Method.Iterations = Method.Iterations + 1; % Iteration count

    %% BEFD Type: iterative method or direct inversion of the system
    % IF one has chosen to use the iterative method
    if (Method.Solver_FD == 0)
    %% Computing ground state using the BEFD method with iterative method
    [FDPhi, flag, relres, iter, resvec] = Local_Krylov_CNFD_solution1d(FDPhi, Method, FDGeometry1D, FDPhysics1D, FDOperators1D); % Computation of the ground state using the iterative BEFD method on a single step of time
    Outputs.iteration_vec(Method.Iterations) = size(resvec,1); % Storing the total number of BICGSTAB iterations at the iteration-th of the Gradient Flow

    % ELSEIF one has chosen to use the direct inversion
    elseif (Method.Solver_FD == 1)
    %% Computing ground state using the BEFD method with direct inversion
    Phi = Local_Direct_CNFD_solution1d(Phi, Method, FDGeometry1D, FDPhysics1D, FDOperators1D); % Computation of the ground state using the direct BEFD method on a single step of time
    end
 
    %% Normalization of the ground states
    Global_L2norm = 0;
    % FOR each component
    for n = 1:Method.Ncomponents
        Global_L2norm = Global_L2norm + L2_norm1d(FDPhi{n},FDGeometry1D)^2; % Computing the L2 norm of each wave function
    end
    % FOR each component
    for n = 1:Method.Ncomponents
        FDPhi{n} = FDPhi{n}/sqrt(Global_L2norm); % Normalization of each wave function
        Method.LocalEvol(n) = max(abs(FDPhi{n}-FDPhi_tmp{n})); % Computing the local evolution of each wave function
    end
    
    %% Updating CPUtime and global evolution
    Method.GlobalEvol = max(Method.LocalEvol); % Computing the global evolution
    Method.Cputime_temp = cputime - Method.Cputime_temp; % Storing the CPUtime relatively to the begining of the "while" loop
    Method.Cputime = Method.Cputime + Method.Cputime_temp; % Updating the CPUtime variable
    
    %% Computation of outputs
    % IF one wants to either compute outputs or to print informations during 
    % the computation
    if (Method.Output) && (mod(Method.Iterations,Outputs.Evo_outputs) ~= 0)
        Outputs = FDOutputs_Var1d(FDPhi, Outputs, Method, FDGeometry1D, FDPhysics1D, FDOperators1D); % Computing the output variables
    end
    %% Printing informations and drawing ground states
    % IF the number of iterations has been reached
    if (mod(Method.Iterations,Print.Evo) == 0)
        % IF one has chosen to print informations
        if (Method.Output) && (Print.Print == 1)
            Print_Info1d(Outputs, Method) % Printing informations
        elseif (~Method.Output) && (Print.Print == 1)
            Outputstmp = FDOutputs_Var1d(FDPhi, Outputs, Method, FDGeometry1D, FDPhysics1D, FDOperators1D); % Computing the output variables
            Print_Info1d(Outputstmp, Method) % Printing informations
        end
        % IF one has chosen to draw the ground states
        if (Print.Draw == 1)
            Draw_solution1d(FDPhi,Method,FDGeometry1D,Figure) % Drawing the wave functions' modulus square and angle
        end;
    end
end

%% Printing the informations of the final ground states
Outputs = FDOutputs_Var1d(FDPhi, Outputs, Method, FDGeometry1D, FDPhysics1D, FDOperators1D); % Computing the output variables
Print_Info1d(Outputs, Method) % Printing informations

%% Finalization of the ground states
% FOR each component
for n = 1:Method.Ncomponents
Phi_final{n} = zeros(Geometry1D.Nx); % Initializing the variable for the storing of the final ground states
Phi_final{n}(2:FDGeometry1D.Nx+1) = FDPhi{n}; % Storing of the ground states solutions
end