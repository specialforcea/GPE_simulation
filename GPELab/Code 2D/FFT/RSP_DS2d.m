%% Relaxation SPectral for the Dynamic System (RSP-DS) : Semi-implicit FFT scheme for the computation of solutions of the dynamic Gross Pitaevskii equation
%% INPUTS:
%%          Phi_0: Initial wave functions (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          Geometry2D: Structure containing variables concerning the geometry of the problem in 2D (structure) (see Geometry2D_Var2d.m)
%%          Physics2D: Structure containing variables concerning the physics of the problem in 2D (structure) (see Physics2D_Var2d.m)
%%          Outputs: Different outputs computed during the computation of the problem (structure) (see OutputsINI_Var2D.m and FFTOutputs_Var2d.m)
%%          Print: Structure containing variables concerning the printing and drawing of informations during the computations (structure) (see Print_Var2d.m)
%%          Figure: Structure containing variables concerning the figures (structure)
%% OUTPUTS:
%%          Phi_final: Dynamic wave functions computated with the CNSP method (cell array)
%%          Outputs: Different outputs computed during the computation of the problem (structure) (see OutputsINI_Var2d.m and FFTOutputs_Var2d.m)
%% FUNCTIONS USED:
%%          FFTGeometry2D_Var2d: To compute the 2d geometry used for the FFT (line 29)
%%          FFTPhysics2D_Var2d: To compute the 2d physics used for the FFT (line 30)
%%          FFTOperators2D_Var2d: To compute the 2d operators used for the FFT (line 38)
%%          Local_Full_CNSP_solution2d: To compute the CNGF on a single time step using the CNSP method (line 53)
%%          FFTOutputs_Var2d: To compute output variables (line 64,68 and 90)
%%          Print_Info2d: To print informations concerning the wave functions (line 71 and 91)
%%          Draw_solution2d: To draw the ground states (line 75)

function [Phi_final, Outputs] = RSP_DS2d(Phi_0, Method, Geometry2D, Physics2D, Outputs, Print, Figure)

%% CPUtime initilization
Method.Cputime = 0; % Initialization of the CPUtime variable
Method.Cputime_temp = cputime; % Initialization of the relative CPUtime variable

%% Geometry, operators, physics and ground states initialization for FFT
FFTGeometry2D = FFTGeometry2D_Var2d(Geometry2D); % Changing the geometry for the FFT
FFTOperators2D = FFTOperators2D_Var2d(FFTGeometry2D); % Computing the derivative FFT operators
% Changing the solution for the FFT
% FOR each component
for n = 1:Method.Ncomponents
    FFTPhi{n} = Phi_0{n}(1:FFTGeometry2D.Ny,1:FFTGeometry2D.Nx); % Removing boundaries for the FFT
end
FFTPhysics2D = FFTPhysics2D_Var2d(FFTPhi, Method, Physics2D, FFTGeometry2D, FFTOperators2D); % Changing the physics for the FFT

%% Initialization of the relaxation variable
Psi = cell(Method.Ncomponents);
FFTPsi = cell(Method.Ncomponents);
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        Psi{n,m} = 0; % Initializing the relaxation variable corresponding to the nonlinearity
    end
end
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        FFTPsi{n,m} = 0; % Initializing the relaxation variable corresponding to the non local nonlinearity
    end
end
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = FFTPhysics2D.Nonlinearity_function_Index{n}
        Psi{n,m} = FFTPhysics2D.Nonlinearity_function{n,m}(FFTPhi,FFTGeometry2D.X,FFTGeometry2D.Y); % Initializing the relaxation variable corresponding to the nonlinearity
    end
end
for n = 1:Method.Ncomponents
    % FOR each component
    for m = FFTPhysics2D.FFTNonlinearity_function_Index{n}
        FFTPsi{n,m} = FFTPhysics2D.FFTNonlinearity_function{n,m}(FFTPhi,FFTGeometry2D.X,FFTGeometry2D.Y,-1i*FFTOperators2D.Gx,-1i*FFTOperators2D.Gy); % Initializing the relaxation variable corresponding to the non local nonlinearity
    end
end

%% CPUtime initilization bis 
Method.Cputime_temp = cputime - Method.Cputime_temp; % Storing the CPUtime relatively to the begining of the program
Method.Cputime  = Method.Cputime + Method.Cputime_temp;  % Updating the CPUtime variable

%% Computation of dynamic of the GPE via RSP
% Stopping criterions: stopping time 
while (Method.Iterations*Method.Deltat<Method.Stop_time)
    %% Updating variables
    FFTPhi_tmp = FFTPhi; % Storing a temporary variable of the ground states to compute local evolution
    Method.Cputime_temp = cputime; % Reinitialization of the relative CPUtime variable
    Method.Iterations = Method.Iterations + 1; % Iteration count
    %% RSP Type: full RSP or partial RSP
    % IF one has chosen to use the full BESP method
    if (Method.Solver_BESP == 0)
    %% Computing ground state using the full RSP method  
    [FFTPhi, Psi, FFTPsi, flag, relres, iter, resvec] = Local_Full_RSP_solution2d(FFTPhi, Psi, FFTPsi, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D); % Computation of the ground state using the CNSP-CNFG method on a single step of time
    Outputs.iteration_vec(Method.Iterations) = size(resvec,1); % Storing the total number of BICGSTAB iterations at the iteration-th iteration of the Gradient Flow
    % ELSEIF one has chosen to use the partial BESP method
    elseif (Method.Solver_BESP == 1)
    [FFTPhi, Psi, FFTPsi,iter] = Local_Partial_RSP_solution2d(FFTPhi, Psi, FFTPsi, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D); % Computation of the ground state using the CNSP-CNFG method on a single step of time
    Outputs.iteration_vec(Method.Iterations) = iter; % Storing the total number of BICGSTAB iterations at the iteration-th iteration of the Gradient Flow
    end
    %% Updating CPUtime
    Method.Cputime_temp = cputime - Method.Cputime_temp; % Storing the CPUtime relatively to the begining of the "while" loop
    Method.Cputime = Method.Cputime + Method.Cputime_temp; % Updating the CPUtime variable

    %% Computation of outputs
    % IF one wants to either compute outputs or to print informations during 
    % the computation
    if (Method.Output) && (mod(Method.Iterations,Outputs.Evo_outputs) == 0)
        Outputs = FFTOutputs_Var2d(FFTPhi, Outputs, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D); % Computing the output variables
    end
    %% Printing informations and drawing ground states
    % IF the number of iterations has been reached
    if (mod(Method.Iterations,Print.Evo) == 0)
        % IF one has chosen to print informations
        if (Method.Output) && (Print.Print == 1)
            Print_Info2d(Outputs, Method) % Printing informations
        elseif (~Method.Output) && (Print.Print == 1)
            Outputstmp = FFTOutputs_Var2d(FFTPhi, Outputs, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D); % Computing the output variables
            Print_Info2d(Outputstmp, Method) % Printing informations
        end
        % IF one has chosen to draw the ground states
        if (Print.Draw == 1)
            Draw_solution2d(FFTPhi,Method,FFTGeometry2D,Figure) % Drawing the wave functions' modulus square and angle
        end;
    end
end

%% Printing the informations of the final ground states
Outputs = FFTOutputs_Var2d(FFTPhi, Outputs, Method, FFTGeometry2D, FFTPhysics2D, FFTOperators2D); % Computing the output variables
Print_Info2d(Outputs, Method) % Printing informations

%% Finalization of the ground states
% FOR each component
for n = 1:Method.Ncomponents
Phi_final{n} = zeros(Geometry2D.Ny,Geometry2D.Nx); % Initializing the variable for the storing of the final ground states
Phi_final{n}(1:FFTGeometry2D.Ny,1:FFTGeometry2D.Nx) = FFTPhi{n}; % Storing of the ground states solutions
Phi_final{n}(FFTGeometry2D.Ny+1,1:FFTGeometry2D.Nx) = Phi_final{n}(1,1:FFTGeometry2D.Nx); % Setting the boundary of the ground states solutions
Phi_final{n}(1:FFTGeometry2D.Ny,FFTGeometry2D.Nx+1) = Phi_final{n}(1:FFTGeometry2D.Ny,1); % Setting the boundary of the ground states solutions
Phi_final{n}(FFTGeometry2D.Ny+1,FFTGeometry2D.Nx+1) = Phi_final{n}(1,1); % Setting the boundary of the ground states solutions
end