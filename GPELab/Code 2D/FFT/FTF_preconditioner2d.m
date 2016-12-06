%% Applying the full Thomas-Fermi preconditioner
%% INPUTS:
%%          Phi_in: Initial component functions (vector)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          FFTGeometry2D: Structure containing variables concerning the geometry of the problem in 2D in the FFT context (structure) (see FFTGeometry2D_Var2d.m)
%%          FFTPhysics2D: Structure containing variables concerning the physics of the problem in 2D in the FFT context (structure) (see FFTPhysics2D_Var2d.m)
%%          FFTOperators2D: Structure containing the derivative FFT operators (structure) (see FFTOperators2D_Var2d.m)
%% OUTPUT:
%%          Phi_out: Component functions with the operators applied (vector)

function [Phi_out] = FTF_preconditioner2d(Phi_in, Method, FFTGeometry2D, FFTPhysics2D)
%% Initialization of variables
Phi_in = reshape(Phi_in,FFTGeometry2D.Ny,Method.Ncomponents*FFTGeometry2D.Nx); % Reshaping vector as a matrix
Phi_out = Phi_in; % Initializing the variable for the component functions with the preconditioner applied

%% Applying the Thomas-Fermi preconditionner
% FOR each component
for n = 1:Method.Ncomponents
    Phi_tmp{n} = Phi_in(:,(1+(n-1)*FFTGeometry2D.Nx):(n*FFTGeometry2D.Nx)); % Exctraction the wave function of a component
end
% FOR each component
for n = 1:Method.Ncomponents
    Phi = zeros(FFTGeometry2D.Ny,FFTGeometry2D.Nx); % Initialization of the wave function
    % FOR each component
    for m = 1:Method.Ncomponents
        Phi = Phi +  FFTPhysics2D.FPThomasFermi{n,m}.*Phi_tmp{m}; % Applying the Thomas-Fermi preconditoner
    end
Phi_out(:,(1+(n-1)*FFTGeometry2D.Nx):(n*FFTGeometry2D.Nx)) = Phi; % Storing the wave function of a component with the Thomas-Fermi preconditoner applied
end

%% Reshapping as a vector the output
Phi_out = reshape(Phi_out,Method.Ncomponents*FFTGeometry2D.N2,1); % Reshapping the wave functions as a vector