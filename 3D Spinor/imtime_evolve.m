function phi_up = imtime_evolve(phi, potential, Deltat, X,Y, Nx,Ny, deltax,deltay,deltafx,deltafy,c0,paritx,parity,dispersion)

mat_pot = potential(X,Y,Z);
mat_nonlin_pot = c0*phi.*conj(phi);

diag_pot = -Deltat*(mat_pot + mat_nonlin_pot);

exp_pot = exp(diag_pot);

phi_up = exp_pot.*phi;
%phi_up = phi_up/norm1d(phi_up, Nx, DeltaX);

disper = dispersion(:,:,:,1);
evol = exp(disper.*Deltat);

fourier_phi = fourier_transform3(phi_up,paritx(:,:,:,1),parity(:,:,:,1),paritz(:,:,:,1),deltax,deltay,deltaz);
%fourier_phi = fourier_phi./norm1d(fourier_phi,Nx,DeltaX);
fourier_phi_evo = evol.*fourier_phi;

phi_up = inverse_ft3(fourier_phi_evo,paritx(:,:,:,1),parity(:,:,:,1),parity(:,:,:,1),deltafx,deltafy,deltafx,Nx,Ny,Nz);

phi_up = phi_up/norm3d(phi_up, Nx, Ny, deltax,deltay);
end


