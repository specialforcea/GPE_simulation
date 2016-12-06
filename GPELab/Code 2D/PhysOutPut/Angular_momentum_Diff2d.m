%% Computation of the angular momentum using the Fourier transform
%% INPUTS:
%%          FDPhi_vect: Wave functions in form of a vector (vector)
%%          FDGeometry2D: Structure containing variables concerning the geometry of the problem in 2D in the FD context (structure) (see FDGeometry2D_Var2d.m)
%%          FDOperators2D: Structure containing the derivative FFT operators (structure) (see FDOperators2D_Var2d.m)
%% OUTPUT:
%%          Angular_momentum: Angular momentum of the wave function (double)

function [Ang] = Angular_momentum_Diff2d(FDPhi, Method, FDGeometry2D, FDOperators2D)
%% Reshaping the wave functions variables as a single vector
FDPhi_vect = zeros(Method.Ncomponents*FDGeometry2D.N2,1); % Initializing the vector to store each component's wave function
% FOR each component
for n = 1:Method.Ncomponents
    FDPhi_vect((1+(n-1)*FDGeometry2D.N2):(n*FDGeometry2D.N2)) = reshape(FDPhi{n},FDGeometry2D.N2,1); % Storing the component's wave function in the vector
end

%% Computation of the momentum and integration over space
Local_Ang =(FDOperators2D.L_z*FDPhi_vect).*conj(FDPhi_vect); %Computing the local angular momentum on the grid
% FOR each component
for n = 1:Method.Ncomponents
Ang{n} = FDGeometry2D.dx*FDGeometry2D.dy*sum(sum(Local_Ang((1+(n-1)*FDGeometry2D.N2):(n*FDGeometry2D.N2)))); % Computing the angular momentum by integrating over space
end