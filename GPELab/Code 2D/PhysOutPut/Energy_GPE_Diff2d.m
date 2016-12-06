%% Computation of the energy of the GPE using the FD operators
%% INPUTS:
%%          Phi: Wave functions (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          FDGeometry2D: Structure containing variables concerning the geometry of the problem in 2D in the FD context (structure) (see FDGeometry2D_Var2d.m)
%%          FDPhysics2D: Structure containing variables concerning the physics of the problem in 2D in the FD context (structure) (see FDPhysics2D_Var2d.m)
%%          FDOperators2D: Structure containing the derivative FD operators (structure) (see FDOperators2D_Var2d.m)
%% OUTPUT:
%%          Energy: Energy of the wave functions (cell array)
%% FUNCTIONS USED:
%%          L2_norm2d: To integrate the local energy (line 54)

function Energy = Energy_GPE_Diff2d(FDPhi, Method, FDGeometry2D, FDPhysics2D, FDOperators2D)
%% Reshaping the wave functions variables as a single vector
FDPhi_vect = zeros(Method.Ncomponents*FDGeometry2D.N2,1); % Initializing the vector to store each component's wave function
% FOR each component
for n = 1:Method.Ncomponents
    FDPhi_vect((1+(n-1)*FDGeometry2D.N2):(n*FDGeometry2D.N2)) = reshape(FDPhi{n},FDGeometry2D.N2,1); % Storing the component's wave function in the vector
end

%% Computing the explicit nonlinearity
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        % IF the energy's nonlinearity is defined
        if (isempty(FDPhysics2D.Nonlinearity_energy_function{n,m}) == 0)
            Nonlinearity_tmp = FDPhysics2D.Nonlinearity_energy_function{n,m}(FDPhi,FDGeometry2D.X,FDGeometry2D.Y); % Computing the nonlinearity
        % ELSE if the energy's nonlinearity is not defined
        else
            Nonlinearity_tmp = FDPhysics2D.Nonlinearity_function{n,m}(FDPhi,FDGeometry2D.X,FDGeometry2D.Y); % Computing the nonlinearity 
        end
        % IF the nonlinearity is a scalar
        if (isscalar(Nonlinearity_tmp))
            NL_scal = Nonlinearity_tmp; % Storing the scalar value of the nonlinearity
            Nonlinearity_tmp = NL_scal*ones(FDGeometry2D.Ny,FDGeometry2D.Nx); % Reconstructing for the nonlinearity
        end
        % IF this is the coupled nonlinearity with the first component
        if (m == 1)
            Nonlinearity_matrix_tmp = spdiags(reshape(Nonlinearity_tmp,FDGeometry2D.N2,1),0,FDGeometry2D.N2,FDGeometry2D.N2); % Initialiazing the temporary local nonlinearity matrix
        % ELSE if this is the coupled nonlinearity with other components
        else
            Nonlinearity_matrix_tmp = horzcat(Nonlinearity_matrix_tmp, spdiags(reshape(Nonlinearity_tmp,FDGeometry2D.N2,1),0,FDGeometry2D.N2,FDGeometry2D.N2)); % Assembling the local nonlinearity matrix for the 'n' component in the 2D geometry for the FD
        end
    end
    % IF this is the first component
    if (n == 1)
        FDOperators2D.Nonlinearity_matrix = Nonlinearity_matrix_tmp;
    % ELSE if this is not the first component
    else
        FDOperators2D.Nonlinearity_matrix = vertcat(FDOperators2D.Nonlinearity_matrix,Nonlinearity_matrix_tmp); % Adding the part of the local nonlinearity matrix corresponding to the 'n' component to the global matrix
    end
    % FOR each component
    for m = 1:Method.Ncomponents
        TimePotential_tmp = FDPhysics2D.TimePotential_function{n,m}(Method.Iterations*Method.Deltat,FDGeometry2D.X,FDGeometry2D.Y); % Computing the time-dependent potential
        % IF the time-dependent potential is a scalar
        if (isscalar(TimePotential_tmp))
            TP_scal = TimePotential_tmp; % Storing the scalar value of the time-dependent potential
            TimePotential_tmp = TP_scal*ones(FDGeometry2D.Ny,FDGeometry2D.Nx); % Reconstructing for the time-dependent potential
        end
        % IF this is the coupled time-dependent potential with the first component
        if (m == 1)
            TimePotential_matrix_tmp = spdiags(reshape(TimePotential_tmp,FDGeometry2D.N2,1),0,FDGeometry2D.N2,FDGeometry2D.N2); % Initialiazing the temporary local time-dependent potential matrix
        % ELSE if this is the coupled time-dependent potential with other components
        else
            TimePotential_matrix_tmp = horzcat(TimePotential_matrix_tmp, spdiags(reshape(TimePotential_tmp,FDGeometry2D.N2,1),0,FDGeometry2D.N2,FDGeometry2D.N2)); % Assembling the local time-dependent potential matrix for the 'n' component in the 2D geometry for the FD
        end
    end
    % IF this is the first component
    if (n == 1)
        FDOperators2D.TimePotential_matrix = TimePotential_matrix_tmp;
    % ELSE if this is not the first component
    else
        FDOperators2D.TimePotential_matrix = vertcat(FDOperators2D.TimePotential_matrix,TimePotential_matrix_tmp); % Adding the part of the local time-dependent potential matrix corresponding to the 'n' component to the global matrix
    end
end

%% Computing energy function
Local_energy = real(conj(FDPhi_vect).*((FDOperators2D.Linear_Operators + FDPhysics2D.Beta*FDOperators2D.Nonlinearity_matrix + FDOperators2D.TimePotential_matrix)*FDPhi_vect)); %Computing the local energy of the GPE

for n = 1:Method.Ncomponents
    %% Integrating over space
    Energy{n} = FDGeometry2D.dx*FDGeometry2D.dy*sum(Local_energy((1+(n-1)*FDGeometry2D.N2):(n*FDGeometry2D.N2))); %Integration of the local energy
end