%% Creation of the 2D operators structure for the FD
%% INPUT:
%%          FDGeometry2D: Structure containing variables concerning the 2D geometry for the FD (structure) (see FDGeometry2D_Var2d.m)
%% OUTPUT:
%%          FDOperators2D: Structure containing the derivative FD operators (structure)
%% FUNCTIONS USED:
%%          Five_Points_Laplacian2d: To compute the second derivative FD operators (line 14)
%%          Two_Points_Rotation2d: To compute the rotational FD operators (line 15)
%%          Two_Points_Gradientx2d: To compute the first derivative in the x direction FD operators (line 16)
%%          Two_Points_Gradienty2d: To compute the first derivative in the y direction FD operators (line 17)

function FDOperators2D = FDOperators2D_Var2d(Method, FDGeometry2D, FDPhysics2D)
%% Computing derivative and identity operators for the FD
Laplacian = Ten_Points_Laplacian2d(FDGeometry2D); %Computing the Laplacian operator discretized with the FD
L_z = Two_Points_Rotation2d(FDGeometry2D); % Computing the rotational operator discretized with the FD
Gx = Four_Points_Gradientx2d(FDGeometry2D); % Computing the gradient operator discretized with the FD in the x direction
Gy = Four_Points_Gradienty2d(FDGeometry2D); % Computing the gradient operator discretized with the FD in the y direction
FDOperators2D.Id = speye(Method.Ncomponents*FDGeometry2D.N2); % Computing the identity operator

%% Computing the multicomponent discretized laplacian and rotational operators
FDOperators2D.L_z = []; % Initializing the discretized rotational operator
FDOperators2D.Laplacian = []; % Initializing the discretized laplacian operator
% FOR each component
for n = 1:Method.Ncomponents
   FDOperators2D.L_z = blkdiag(FDOperators2D.L_z , L_z); % Adding a local discretized rotational operator
   FDOperators2D.Laplacian = blkdiag(FDOperators2D.Laplacian , Laplacian); % Adding a local discretized laplacian operator
end

%% Computing the potential and gradients operators in the 2D geometry for the FD
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        % IF the poetntial is a scalar
        if (isscalar(FDPhysics2D.Potential{n,m}))
            Potential_scal = FDPhysics2D.Potential{n,m}; % Storing the scalar value of the gradient
            FDPhysics2D.Potential{n,m} = Potential_scal*ones(FDGeometry2D.Ny,FDGeometry2D.Nx); % Reconstructing for the gradient
        % ELSEIF the potential is empty
        elseif (isempty(FDPhysics2D.Potential{n,m}))
            FDPhysics2D.Potential{n,m} = zeros(FDGeometry2D.Ny,FDGeometry2D.Nx); % Reconstructing for the gradient
        end
        % IF this is the coupled potential with the first component
        if (m == 1)
            Potential_matrix_tmp = spdiags(reshape(FDPhysics2D.Potential{n,m},FDGeometry2D.N2,1),0,FDGeometry2D.N2,FDGeometry2D.N2); % Initialiazing the temporary local potential matrix
        % ELSE if this is the coupled potential with other components
        else
            Potential_matrix_tmp = horzcat(Potential_matrix_tmp, spdiags(reshape(FDPhysics2D.Potential{n,m},FDGeometry2D.N2,1),0,FDGeometry2D.N2,FDGeometry2D.N2)); % Assembling the local potential matrix for the 'n' component in the 2D geometry for the FD
        end
    end
    % IF this is the first component
    if (n == 1)
        FDOperators2D.Potential_matrix = Potential_matrix_tmp;
    % ELSE if this is not the first component
    else
        FDOperators2D.Potential_matrix = vertcat(FDOperators2D.Potential_matrix,Potential_matrix_tmp); % Adding the part of the local potential matrix corresponding to the 'n' component to the global matrix
    end
end

% Computing the gradient in the x direction functions
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        % IF the gradient is a scalar
        if (isscalar(FDPhysics2D.Gradientx{n,m}))
            Grad_scal = FDPhysics2D.Gradientx{n,m}; % Storing the scalar value of the gradient
            FDPhysics2D.Gradientx{n,m} = Grad_scal*ones(FDGeometry2D.Ny,FDGeometry2D.Nx); % Reconstructing for the gradient
        % ELSEIF the gradient is empty
        elseif (isempty(FDPhysics2D.Gradientx{n,m}))
            FDPhysics2D.Gradientx{n,m} = zeros(FDGeometry2D.Ny,FDGeometry2D.Nx); % Reconstructing for the gradient
        end
        % IF this is the coupled gradient in the x direction with the first component
        if (m == 1)
            Gradientx_matrix_tmp = spdiags(reshape(FDPhysics2D.Gradientx{n,m},FDGeometry2D.N2,1),0,FDGeometry2D.N2,FDGeometry2D.N2)*Gx; % Initialiazing the temporary local gradient in the x direction matrix
        % ELSE if this is the coupled gradient in the x direction with other components
        else
            Gradientx_matrix_tmp = horzcat(Gradientx_matrix_tmp, spdiags(reshape(FDPhysics2D.Gradientx{n,m},FDGeometry2D.N2,1),0,FDGeometry2D.N2,FDGeometry2D.N2)*Gx); % Assembling the local potential gradient in the x direction for the 'n' component in the 2D geometry for the FD
        end
    end
    % IF this is the first component
    if (n == 1)
        FDOperators2D.Gradientx_matrix = Gradientx_matrix_tmp;
    % ELSE if this is not the first component
    else
        FDOperators2D.Gradientx_matrix = vertcat(FDOperators2D.Gradientx_matrix,Gradientx_matrix_tmp); % Adding the part of the local gradient in the x direction matrix corresponding to the 'n' component to the global matrix
    end
end

% Computing the gradient in the y direction functions
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        % IF the gradient is a scalar
        if (isscalar(FDPhysics2D.Gradienty{n,m}))
            Grad_scal = FDPhysics2D.Gradienty{n,m}; % Storing the scalar value of the gradient
            FDPhysics2D.Gradienty{n,m} = Grad_scal*ones(FDGeometry2D.Ny,FDGeometry2D.Nx); % Reconstructing for the gradient
        % ELSEIF the gradient is empty
        elseif (isempty(FDPhysics2D.Gradienty{n,m}))
            FDPhysics2D.Gradienty{n,m} = zeros(FDGeometry2D.Ny,FDGeometry2D.Nx); % Reconstructing for the gradient
        end
        % IF this is the coupled potential with the first component
        if (m == 1)
            Gradienty_matrix_tmp = spdiags(reshape(FDPhysics2D.Gradienty{n,m},FDGeometry2D.N2,1),0,FDGeometry2D.N2,FDGeometry2D.N2)*Gy; % Initialiazing the temporary local gradient in the y direction matrix
        % ELSE if this is the coupled potential with other components
        else
            Gradienty_matrix_tmp = horzcat(Gradienty_matrix_tmp, spdiags(reshape(FDPhysics2D.Gradienty{n,m},FDGeometry2D.N2,1),0,FDGeometry2D.N2,FDGeometry2D.N2)*Gy); % Assembling the local potential gradient in the y direction for the 'n' component in the 2D geometry for the FD
        end
    end
    % IF this is the first component
    if (n == 1)
        FDOperators2D.Gradienty_matrix = Gradienty_matrix_tmp;
    % ELSE if this is not the first component
    else
        FDOperators2D.Gradienty_matrix = vertcat(FDOperators2D.Gradienty_matrix,Gradienty_matrix_tmp); % Adding the part of the local gradient in the y direction matrix corresponding to the 'n' component to the global matrix
    end
end

FDOperators2D.Linear_Operators = FDOperators2D.Potential_matrix + FDOperators2D.Gradientx_matrix + FDOperators2D.Gradienty_matrix - FDPhysics2D.Delta*FDOperators2D.Laplacian;