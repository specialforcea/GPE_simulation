%% Applying the Thomas-Fermi preconditioner
%% INPUTS:
%%          Phi_in: Initial component functions (vector)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var3d.m)
%%          FFTGeometry3D: Structure containing variables concerning the geometry of the problem in 3D in the FFT context (structure) (see FFTGeometry3D_Var3d.m)
%%          FFTPhysics3D: Structure containing variables concerning the physics of the problem in 3D in the FFT context (structure) (see FFTPhysics3D_Var3d.m)
%%          FFTOperators3D: Structure containing the derivative FFT operators (structure) (see FFTOperators3D_Var3d.m)
%% OUTPUT:
%%          Phi_out: Component functions with the operators applied (vector)

function [Phi_out] = TF_preconditioner3d(Phi_in, Method, FFTGeometry3D, FFTPhysics3D)
%% Initialization of variables
Phi_in = reshape(Phi_in,FFTGeometry3D.Ny,FFTGeometry3D.Nx,Method.Ncomponents*FFTGeometry3D.Nz); % Reshaping vector as a matrix
Phi_out = Phi_in; % Initializing the variable for the component functions with the preconditionner applied

%% Applying the Laplace preconditioner
% FOR each component
for n = 1:Method.Ncomponents
Phi = Phi_in(:,:,(1+(n-1)*FFTGeometry3D.Nz):(n*FFTGeometry3D.Nz)); % Exctraction the wave function of a component
Phi = 1./(1/Method.Deltat + FFTPhysics3D.Potential{n,n} + FFTPhysics3D.TimePotential{n,n} + FFTPhysics3D.Beta*FFTPhysics3D.Nonlinearity{n,n} + FFTPhysics3D.Beta*FFTPhysics3D.FFTNonlinearity{n,n}).*Phi; % Applying the Laplace preconditoner
Phi_out(:,:,(1+(n-1)*FFTGeometry3D.Nz):(n*FFTGeometry3D.Nz)) = Phi; % Storing the wave function of a component with the Laplace preconditoner applied
end

%% Reshapping as a vector the output
Phi_out = reshape(Phi_out,Method.Ncomponents*FFTGeometry3D.N3,1); % Reshapping the wave functions as a vector