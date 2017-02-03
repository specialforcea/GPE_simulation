function [phi_0] = Thomas_fermi2D(Beta,X,TF_radius,potential,Nx,DeltaX)
phi_0 = real(sqrt(((Beta/pi)^(1/2) - potential)/Beta)); % Computing the Thomas-Fermi approximation
%phi_0 = phi_0./norm1d(phi_0,Nx,DeltaX);
%phi_0 = real(sqrt(pi/(4*Beta).*(TF_radius^2 - X.^2).^2));%integrate over 2d to get 1d wavefunction.

cutX = abs(X) < TF_radius;
phi_0 = phi_0.*cutX;

end