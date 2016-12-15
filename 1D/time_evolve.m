function phi_up = time_evolve(phi, potential, Deltat, X, Beta, Nx, deltax,deltaf,L)
mat_pot = potential(X);
mat_nonlin_pot = 0;%Beta*phi.*conj(phi);

diag_pot = -Deltat*(mat_pot + mat_nonlin_pot);

exp_pot = exp(diag_pot);

phi_up = exp_pot.*phi;
%phi_up = phi_up/norm1d(phi_up, Nx, DeltaX);
lin_n =  (0:1:Nx-1) - (Nx-1)/2;

dispersion = -(2*pi/L*lin_n).^2 *Deltat/2;
evol = exp(dispersion);

fourier_phi = fourier_transform(phi_up,Nx,deltax);
%fourier_phi = fourier_phi./norm1d(fourier_phi,Nx,DeltaX);
fourier_phi_evo = evol.*fourier_phi;

phi_up = inverse_ft(fourier_phi_evo,Nx,deltaf);

phi_up = phi_up/norm1d(phi_up, Nx, deltax);
end


