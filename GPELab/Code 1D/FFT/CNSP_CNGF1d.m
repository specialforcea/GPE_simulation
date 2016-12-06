%% Crank Nichoslon SPectral for the Continuous Normalized Gradient Flow (CNSP-CNFG) : Semi-implicit FFT scheme for the computation of grounds states of the Gross Pitaevskii equation
%% INPUTS:
%%          Phi_0: Initial wave functions (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var1d.m)
%%          Geometry1D: Structure containing variables concerning the geometry of the problem in 1D (structure) (see Geometry1D_Var1d.m)
%%          Physics1D: Structure containing variables concerning the physics of the problem in 1D (structure) (see Physics1D_Var1d.m)
%%          Outputs: Different outputs computed during the computation of the ground state (structure) (see OutputsINI_Var1D.m and FFTOutputs_Var1d.m)
%%          Print: Structure containing variables concerning the printing and drawing of informations during the computations (structure) (see Print_Var1d.m)
%%          Figure: Structure containing variables concerning the figures (structure)
%% OUTPUTS:
%%          Phi_final: Ground states computated with the CNSP method (cell array)
%%          Outputs: Different outputs computed during the computation of the ground state (structure) (see OutputsINI_Var1d.m and FFTOutputs_Var1d.m)
%% FUNCTIONS USED:
%%          FFTGeometry1D_Var1d: To compute the 1d geometry used for the FFT (line 30)
%%          FFTPhysics1D_Var1d: To compute the 1d physics used for the FFT (line 31)
%%          FFTOperators2D_Var2d: To compute the 1d operators used for the FFT (line 39)
%%          Local_Full_CNSP_solution1d: To compute the CNGF on a single time step using the CNSP method (line 59)
%%          L2_norm1d: To compute the normalized wave function (line 65)
%%          FFTOutputs_Var1d: To compute output variables (line 78,82 and 95)
%%          Print_Info1d: To print informations concerning the wave functions (line 85 and 96)
%%          Draw_solution1d: To draw the ground states (line 89)

function [Phi_final, Outputs] = CNSP_CNGF1d(Phi_0, Method, Geometry1D, Physics1D, Outputs, Print, Figure)

%% CPUtime initilization
Method.Cputime = 0; % Initialization of the CPUtime variable
Method.Cputime_temp = cputime; % Initialization of the relative CPUtime variable

%% Geometry, operators, physics and ground states initialization for FFT
FFTGeometry1D = FFTGeometry1D_Var1d(Geometry1D); % Changing the geometry for the FFT
FFTOperators1D = FFTOperators1D_Var1d(FFTGeometry1D); % Computing the derivative FFT operators
% Changing the solution for the FFT
% FOR each component
for n = 1:Method.Ncomponents
    FFTPhi{n} = Phi_0{n}(1:FFTGeometry1D.Nx); % Removing boundaries for the FFT
end
FFTPhysics1D = FFTPhysics1D_Var1d(FFTPhi, Method, Physics1D, FFTGeometry1D, FFTOperators1D); % Changing the physics for the FFT

%% CPUtime initilization bis 
Method.Cputime_temp = cputime - Method.Cputime_temp; % Storing the CPUtime relatively to the begining of the program
Method.Cputime  = Method.Cputime + Method.Cputime_temp;  % Updating the CPUtime variable

%% Computation of the ground state via CNSP
% Stopping criterions:- the global evolution of the ground states (maximum
%                     of each local evolution) must be lower than the 
%                     stopping criterion Method.Stop multiplied by the
%                     time step
%                     - the number of iterations must not exceed the
%                     maximum number of iterations  
while (Method.EvolutionCriterion > Method.Stop_crit{2}*Method.Deltat)&&(Method.Iterations<Method.Max_iter)
    %% Updating variables
    FFTPhi_tmp = FFTPhi; % Storing a temporary variable of the ground states to compute local evolution
    Method.Cputime_temp = cputime; % Reinitialization of the relative CPUtime variable
    Method.Iterations = Method.Iterations + 1; % Iteration count

    %% Computing ground state using the full BESP method  
        [FFTPhi,flag, relres, iter, resvec] = Local_Full_CNSP_solution1d(FFTPhi, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D); % Computation of the ground state using the CNSP-CNFG method on a single step of time
        Outputs.iteration_vec(Method.Iterations) = size(resvec,1); % Storing the total number of BICGSTAB iterations at the iteration-th iteration of the Gradient Flow

    %% Normalization of the ground states
    if strcmp(Method.Normalization,'Multi')
    Global_L2norm = 0;
    for n = 1:Method.Ncomponents
        Global_L2norm = Global_L2norm + L2_norm1d(FFTPhi{n},FFTGeometry1D)^2; % Computing the norm of each wave function
    end
    for n = 1:Method.Ncomponents
        FFTPhi{n} = FFTPhi{n}/sqrt(Global_L2norm)*sqrt(Method.NParticles(1)); % Normalization of each wave function
        Method.LocalEvol(n) = max(abs(FFTPhi{n}-FFTPhi_tmp{n})); % Computing the local evolution of each wave function
    end
    elseif strcmp(Method.Normalization,'Single')
    for n = 1:Method.Ncomponents
        FFTPhi{n} = FFTPhi{n}/L2_norm1d(FFTPhi{n},FFTGeometry1D)*sqrt(Method.NParticles(n)); % Normalization of each wave function
        Method.LocalEvol(n) = max(abs(FFTPhi{n}-FFTPhi_tmp{n})); % Computing the local evolution of each wave function
    end
    end
    
    %% Updating CPUtime and evolution criterions
    if strcmp(Method.Stop_crit{1},'MaxNorm')
        Method.EvolutionCriterion = max(Method.LocalEvol); % Computing the global evolution
    elseif strcmp(Method.Stop_crit{1},'Energy')
        Energy = Energy_GPE_Fourier1d(FFTPhi, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D); % Computing the energy of each wave function
        if (Method.Iterations == 1)
            Energy_tmp = Addc(Energy,1);
        end
        EnergyEvolution = 0; % Initializing the energy evolution variable
        for n = 1:Method.Ncomponents
            EnergyEvolution = EnergyEvolution + abs(Energy{n}-Energy_tmp{n});
        end
        Energy_tmp = Energy;
        Method.EvolutionCriterion = EnergyEvolution;
    end
    Method.Cputime_temp = cputime - Method.Cputime_temp; % Storing the CPUtime relatively to the begining of the "while" loop
    Method.Cputime = Method.Cputime + Method.Cputime_temp; % Updating the CPUtime variable

    %% Computation of outputs
    % IF one wants to either compute outputs or to print informations during 
    % the computation
    if (Method.Output) && (mod(Method.Iterations,Outputs.Evo_outputs) == 0)
        Outputs = FFTOutputs_Var1d(FFTPhi, Outputs, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D); % Computing the output variables
    end
    %% Printing informations and drawing ground states
    % IF the number of iterations has been reached
    if (Method.Output) && (mod(Method.Iterations,Print.Evo) == 0)
        % IF one has chosen to print informations
        if (Method.Output) && (Print.Print == 1)
            Print_Info1d(Outputs, Method) % Printing informations
        elseif (~Method.Output) && (Print.Print == 1)
            Outputstmp = FFTOutputs_Var1d(FFTPhi, Outputs, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D); % Computing the output variables
            Print_Info1d(Outputstmp, Method) % Printing informations
        end
        % IF one has chosen to draw the ground states
        if (Print.Draw == 1)
            Draw_solution1d(FFTPhi,Method,FFTGeometry1D,Figure) % Drawing the wave functions' modulus square and angle
        end;
    end
end

%% Printing the informations of the final ground states
Outputs = FFTOutputs_Var1d(FFTPhi, Outputs, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D); % Computing the output variables
Print_Info1d(Outputs, Method) % Printing informations

%% Finalization of the ground states
% FOR each component
for n = 1:Method.Ncomponents
Phi_final{n} = zeros(Geometry1D.Nx,1); % Initializing the variable for the storing of the final ground states
Phi_final{n}(1:FFTGeometry1D.Nx) = FFTPhi{n}; % Storing of the ground states solutions
Phi_final{n}(FFTGeometry1D.Nx+1) = Phi_final{n}(1); % Setting the periodic boundary of the ground states solutions
end