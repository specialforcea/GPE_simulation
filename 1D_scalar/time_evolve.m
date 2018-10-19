function phi_up = time_evolve(phi, potential, Deltat, X, Nx, deltax,deltaf,L,c0)

density = sq(phi);
tot_density = sum(density);

pot = potential(X);


nonlin_pot = c0.*tot_density;


diag_pot = -Deltat.*(pot + nonlin_pot);

exp_pot = exp(diag_pot);

phi_up = exp_pot.*phi;

if(mod(Nx,2) ~= 0)
    lin_n =  (0:1:Nx-1) - (Nx-1)/2;
else
    N = Nx + 1;
    lin_n = (0:1:N-1) - (N-1)/2;
    lin_n = lin_n(1:Nx);
end
dispersion = -(2*pi/L*lin_n).^2 *Deltat/2;
evol = exp(dispersion);


fourier_phi = fourier_transform(phi_up,Nx,deltax);


fourier_phi_evo = evol.*fourier_phi;

phi_up = inverse_ft(fourier_phi_evo,Nx,deltaf);


phi_up = phi_up/norm1d(phi_up, Nx, deltax);



end
