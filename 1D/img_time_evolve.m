function phi_up = img_time_evolve(phi, potential, Deltat, X, Beta)
mat_pot = potential(X);
mat_nonlin_pot = Beta*phi.*conj(phi);

diag_pot = -Deltat*(mat_pot + mat_nonlin_pot);

exp_pot = exp(diag_pot);

phi_up = exp_pot*phi;

FFTphi = fft(phi_up);


