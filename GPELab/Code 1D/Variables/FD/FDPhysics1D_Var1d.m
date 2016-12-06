%% Creation of the 1D physics structure for the FD
%% INPUTS:
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var1d.m)
%%          Physics1D: Structure containing variables concerning the physics of the problem in 1D (structure) (see Physics1D_Var1d.m)
%%          FDGeometry1D: Structure containing variables concerning the 1D geometry for the FD (structure) (see FDGeometry1D_Var1d.m)
%% OUTPUT:
%%          FDPhysics1D: Structure containing variables concerning the physics of the problem in 1D in the FD context (structure)

function [FDPhysics1D] = FDPhysics1D_Var1d(Method, Physics1D, FDGeometry1D)
%% Initialization
FDPhysics1D = Physics1D; % Copying the 1D physics

%% Computing the potential and gradients functions in the 1D geometry for the FD
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        FDPhysics1D.Potential{n,m} = Physics1D.Potential_function{n,m}(FDGeometry1D.X); % Computing the potential function in the 1D geometry for the FD
    end
end

% Computing the gradient in the x direction functions
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        FDPhysics1D.Gradientx{n,m} = Physics1D.Gradientx_function{n,m}(FDGeometry1D.X); % Computing the gradient in the x direction function in the 1D geometry for the FD
    end
end