%% Creation of the 2D physics structure for the FD
%% INPUTS:
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          Physics2D: Structure containing variables concerning the physics of the problem in 2D (structure) (see Physics2D_Var2d.m)
%%          FDGeometry2D: Structure containing variables concerning the 2D geometry for the FD (structure) (see FDGeometry2D_Var2d.m)
%% OUTPUT:
%%          FDPhysics2D: Structure containing variables concerning the physics of the problem in 2D in the FD context (structure)

function [FDPhysics2D] = FDPhysics2D_Var2d(Method, Physics2D, FDGeometry2D)
%% Initialization
FDPhysics2D = Physics2D; % Copying the 2D physics

%% Computing the potential and gradients functions in the 2D geometry for the FD
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        FDPhysics2D.Potential{n,m} = Physics2D.Potential_function{n,m}(FDGeometry2D.X,FDGeometry2D.Y); % Computing the potential function in the 2D geometry for the FD
    end
end

% Computing the gradient in the x direction functions
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        FDPhysics2D.Gradientx{n,m} = Physics2D.Gradientx_function{n,m}(FDGeometry2D.X,FDGeometry2D.Y); % Computing the gradient in the x direction function in the 2D geometry for the FD
    end
end

% Computing the gradient in the y direction functions
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component
    for m = 1:Method.Ncomponents
        FDPhysics2D.Gradienty{n,m} = Physics2D.Gradienty_function{n,m}(FDGeometry2D.X,FDGeometry2D.Y); % Computing the gradient in the y direction function in the 2D geometry for the FD
    end
end