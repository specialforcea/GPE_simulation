%% Applying the Thomas-Fermi preconditioner
%% INPUTS:
%%          Phi_in: Initial component functions (vector)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var3d.m)
%%          FFTGeometry3D: Structure containing variables concerning the geometry of the problem in 3D in the FFT context (structure) (see FFTGeometry3D_Var3d.m)
%%          FFTPhysics3D: Structure containing variables concerning the physics of the problem in 3D in the FFT context (structure) (see FFTPhysics3D_Var3d.m)
%% OUTPUT:
%%          Phi_out: Component functions with the operators applied (vector)

function [Phi_out] = FTF_preconditioner3d(Phi_in, Method, FFTGeometry3D, FFTPhysics3D)
%% Initialization of variables
Phi_in = reshape(Phi_in,FFTGeometry3D.Ny,FFTGeometry3D.Nx,Method.Ncomponents*FFTGeometry3D.Nz); % Reshaping vector as a matrix
Phi_out = zeros(FFTGeometry3D.Ny,FFTGeometry3D.Nx,Method.Ncomponents*FFTGeometry3D.Nz); % Initializing the variable for the component functions with the preconditionner applied
% FOR each component
for n = 1:Method.Ncomponents
    Phi{n} = Phi_in(:,:,(1+(n-1)*FFTGeometry3D.Nz):(n*FFTGeometry3D.Nz)); % Exctraction the wave function of a component
end

%% Applying the ThomasFermi preconditioner
for n = 1:Method.Ncomponents
    Phi_tmp = zeros(FFTGeometry3D.Ny,FFTGeometry3D.Nx,FFTGeometry3D.Nz);
    for m = 1:Method.Ncomponents
        Phi_tmp = Phi_tmp + FFTPhysics3D.FThomasFermi{n,m}.*Phi{m}; % Applying the Thomas-Fermi preconditioner
    end
    Phi_out(:,:,(1+(n-1)*FFTGeometry3D.Nz):(n*FFTGeometry3D.Nz)) = Phi_tmp; % Storing the wave function with the operators and the preconditionner applied
end

%% Reshapping as a vector the output
Phi_out = reshape(Phi_out,Method.Ncomponents*FFTGeometry3D.N3,1); % Reshapping the wave functions as a vector