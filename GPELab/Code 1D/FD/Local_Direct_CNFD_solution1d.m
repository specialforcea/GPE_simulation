%% Computation of a step in time using a direct inversion with CNFD scheme
%% INPUTS:
%%          Phi: Initial wave functions in the 2D geometry for the FD (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var1d.m)
%%          FDGeometry1D: Structure containing variables concerning the geometry of the problem in 1D in the FD context (structure) (see FDGeometry1D_Var1d.m)
%%          FDPhysics1D: Structure containing variables concerning the physics of the problem in 1D in the FD context (structure) (see FDPhysics1D_Var1d.m)
%%          FDOperators1D: Structure containing the derivative FFT operators (structure) (see FDOperators1D_Var1d.m)
%% OUTPUTS:
%%          Phi: Wave functions computated with the CNFD method on a single step (cell array)

function Phi = Local_Direct_CNFD_solution1d(Phi, Method, FDGeometry1D, FDPhysics1D, FDOperators1D)

%% Computing the explicit nonlinearity and the time-dependent potential
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        Nonlinearity_tmp = FDPhysics1D.Nonlinearity_function{n,m}(Phi,FDGeometry1D.X); % Computing the nonlinearity
        % IF the nonlinearity is a scalar
        if (isscalar(Nonlinearity_tmp))
            NL_scal = Nonlinearity_tmp; % Storing the scalar value of the nonlinearity
            Nonlinearity_tmp = NL_scal*ones(FDGeometry1D.Nx,1); % Reconstructing for the nonlinearity
        end
        % IF this is the coupled nonlinearity with the first component
        if (m == 1)
            Nonlinearity_matrix_tmp = spdiags(reshape(Nonlinearity_tmp,FDGeometry1D.Nx,1),0,FDGeometry1D.Nx,FDGeometry1D.Nx); % Initialiazing the temporary local nonlinearity matrix
        % ELSE if this is the coupled nonlinearity with other components
        else
            Nonlinearity_matrix_tmp = horzcat(Nonlinearity_matrix_tmp, spdiags(reshape(Nonlinearity_tmp,FDGeometry1D.Nx,1),0,FDGeometry1D.Nx,FDGeometry1D.Nx)); % Assembling the local nonlinearity matrix for the 'n' component in the 2D geometry for the FD
        end
    end
    % IF this is the first component
    if (n == 1)
        FDOperators1D.Nonlinearity_matrix = Nonlinearity_matrix_tmp;
    % ELSE if this is not the first component
    else
        FDOperators1D.Nonlinearity_matrix = vertcat(FDOperators1D.Nonlinearity_matrix,Nonlinearity_matrix_tmp); % Adding the part of the local nonlinearity matrix corresponding to the 'n' component to the global matrix
    end
    % FOR each component
    for m = 1:Method.Ncomponents
        TimePotential_tmp = FDPhysics1D.TimePotential_function{n,m}(Method.Iterations*Method.Deltat,FDGeometry1D.X); % Computing the nonlinearity
        % IF the time-dependent potential is a scalar
        if (isscalar(TimePotential_tmp))
            TP_scal = TimePotential_tmp; % Storing the scalar value of the nonlinearity
            TimePotential_tmp = TP_scal*ones(FDGeometry1D.Nx,1); % Reconstructing for the nonlinearity
        end
        % IF this is the coupled time-dependent potential with the first component
        if (m == 1)
            TimePotential_matrix_tmp = spdiags(reshape(TimePotential_tmp,FDGeometry1D.Nx,1),0,FDGeometry1D.Nx,FDGeometry1D.Nx); % Initialiazing the temporary local time-dependent potential matrix
        % ELSE if this is the coupled nonlinearity with other components
        else
            TimePotential_matrix_tmp = horzcat(TimePotential_matrix_tmp, spdiags(reshape(TimePotential_tmp,FDGeometry1D.Nx,1),0,FDGeometry1D.Nx,FDGeometry1D.Nx)); % Assembling the local time-dependent potential matrix for the 'n' component in the 2D geometry for the FD
        end
    end
    % IF this is the first component
    if (n == 1)
        FDOperators1D.TimePotential_matrix = TimePotential_matrix_tmp;
    % ELSE if this is not the first component
    else
        FDOperators1D.TimePotential_matrix = vertcat(FDOperators1D.TimePotential_matrix,TimePotential_matrix_tmp); % Adding the part of the local time-dependent potential matrix corresponding to the 'n' component to the global matrix
    end
end

%% Reshaping the wave functions variables as a single vector for the direct method
Phi_vect = zeros(Method.Ncomponents*FDGeometry1D.Nx,1); % Initializing the vector to store each component's wave function
% FOR each component
for n = 1:Method.Ncomponents
    Phi_vect((1+(n-1)*FDGeometry1D.Nx):(n*FDGeometry1D.Nx)) = reshape(Phi{n},FDGeometry1D.Nx,1); % Storing the component's wave function in the vector
end

%% Computing the wave functions after a single time step using a direct inversion for CNFD
%If the computation is dynamic
if (strcmp(Method.Computation,'Dynamic'))
    Method.Deltat = 1i*Method.Deltat;
end
Phi_vect = ((1/Method.Deltat)*FDOperators1D.Id - FDOperators1D.Linear_Operators/2 - FDPhysics1D.Beta*FDOperators1D.Nonlinearity_matrix/2 - FDOperators1D.TimePotential_matrix/2)*Phi_vect; % Computation of the explicit system
Phi_vect = ((1/Method.Deltat)*FDOperators1D.Id + FDOperators1D.Linear_Operators/2 + FDPhysics1D.Beta*FDOperators1D.Nonlinearity_matrix/2 + FDOperators1D.TimePotential_matrix/2)\Phi_vect; % Direct inversion of the system

%% Reshaping vector output from the iterative method as matrix
% FOR each component
for n = 1:Method.Ncomponents
    Phi{n} = reshape(Phi_vect((1+(n-1)*FDGeometry1D.Nx):(n*FDGeometry1D.Nx)),FDGeometry1D.Nx,1); % Storing the components' wave functions back in a cell array
end