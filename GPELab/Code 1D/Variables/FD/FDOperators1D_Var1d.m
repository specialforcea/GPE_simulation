%% Creation of the 1D operators structure for the FD
%% INPUT:
%%          FDGeometry1D: Structure containing variables concerning the 1D geometry for the FD (structure) (see FDGeometry1D_Var1d.m)
%% OUTPUT:
%%          FDOperators1D: Structure containing the derivative FD operators (structure)
%% FUNCTIONS USED:
%%          FivePoints_Laplacian1d: To compute the second derivative FD operators (line 12)
%%          Two_Points_Gradientx1d: To compute the first derivative in the x direction FD operators (line 13)

function FDOperators1D = FDOperators1D_Var1d(Method, FDGeometry1D, FDPhysics1D)
%% Computing derivative and identity operators for the FD
Laplacian = Five_Points_Laplacian1d(FDGeometry1D); %Computing the Laplacian operator discretized with the FD
Gx = Four_Points_Gradientx1d(FDGeometry1D); % Computing the gradient operator discretized with the FD in the x direction
FDOperators1D.Id = speye(Method.Ncomponents*FDGeometry1D.Nx); % Computing the identity operator

%% Computing the multicomponent discretized laplacian operator
FDOperators1D.Laplacian = []; % Initializing the discretized laplacian operator
% FOR each component
for n = 1:Method.Ncomponents
   FDOperators1D.Laplacian = blkdiag(FDOperators1D.Laplacian , Laplacian); % Adding a local discretized laplacian operator
end

%% Computing the potential and gradients operators in the 1D geometry for the FD
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        % IF the poetntial is a scalar
        if (isscalar(FDPhysics1D.Potential{n,m}))
            Potential_scal = FDPhysics1D.Potential{n,m}; % Storing the scalar value of the potential
            FDPhysics1D.Potential{n,m} = Potential_scal*ones(FDGeometry1D.Nx,1); % Reconstructing for the potential
        % ELSEIF the potential is empty
        elseif (isempty(FDPhysics1D.Potential{n,m}))
            FDPhysics1D.Potential{n,m} = zeros(FDGeometry1D.Nx,1); % Reconstructing for the potential
        end
        % IF this is the coupled potential with the first component
        if (m == 1)
            Potential_matrix_tmp = spdiags(reshape(FDPhysics1D.Potential{n,m},FDGeometry1D.Nx,1),0,FDGeometry1D.Nx,FDGeometry1D.Nx); % Initialiazing the temporary local potential matrix
        % ELSE if this is the coupled potential with other components
        else
            Potential_matrix_tmp = horzcat(Potential_matrix_tmp, spdiags(reshape(FDPhysics1D.Potential{n,m},FDGeometry1D.Nx,1),0,FDGeometry1D.Nx,FDGeometry1D.Nx)); % Assembling the local potential matrix for the 'n' component in the 2D geometry for the FD
        end
    end
    % IF this is the first component
    if (n == 1)
        FDOperators1D.Potential_matrix = Potential_matrix_tmp;
    % ELSE if this is not the first component
    else
        FDOperators1D.Potential_matrix = vertcat(FDOperators1D.Potential_matrix,Potential_matrix_tmp); % Adding the part of the local potential matrix corresponding to the 'n' component to the global matrix
    end
end

% Computing the gradient in the x direction functions
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        % IF the gradient is a scalar
        if (isscalar(FDPhysics1D.Gradientx{n,m}))
            Grad_scal = FDPhysics1D.Gradientx{n,m}; % Storing the scalar value of the gradient
            FDPhysics1D.Gradientx{n,m} = Grad_scal*ones(FDGeometry1D.Nx,1); % Reconstructing for the gradient
        % ELSEIF the gradient is empty
        elseif (isempty(FDPhysics1D.Gradientx{n,m}))
            FDPhysics1D.Gradientx{n,m} = zeros(FDGeometry1D.Nx,1); % Reconstructing for the gradient
        end
        % IF this is the coupled gradient in the x direction with the first component
        if (m == 1)
            Gradientx_matrix_tmp = spdiags(reshape(FDPhysics1D.Gradientx{n,m},FDGeometry1D.Nx,1),0,FDGeometry1D.Nx,FDGeometry1D.Nx)*Gx; % Initialiazing the temporary local gradient in the x direction matrix
        % ELSE if this is the coupled gradient in the x direction with other components
        else
            Gradientx_matrix_tmp = horzcat(Gradientx_matrix_tmp, spdiags(reshape(FDPhysics1D.Gradientx{n,m},FDGeometry1D.Nx,1),0,FDGeometry1D.Nx,FDGeometry1D.Nx)*Gx); % Assembling the local potential gradient in the x direction for the 'n' component in the 2D geometry for the FD
        end
    end
    % IF this is the first component
    if (n == 1)
        FDOperators1D.Gradientx_matrix = Gradientx_matrix_tmp;
    % ELSE if this is not the first component
    else
        FDOperators1D.Gradientx_matrix = vertcat(FDOperators1D.Gradientx_matrix,Gradientx_matrix_tmp); % Adding the part of the local gradient in the x direction matrix corresponding to the 'n' component to the global matrix
    end
end

FDOperators1D.Linear_Operators = FDOperators1D.Potential_matrix + FDOperators1D.Gradientx_matrix - FDPhysics1D.Delta*FDOperators1D.Laplacian;