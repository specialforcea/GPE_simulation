function y = chem_pot(phi,X,Nx,Beta,k_scale,deltax,deltaf,V,L)
lin_n =  (0:1:Nx-1) - (Nx-1)/2;
fourier_phi = fourier_transform(phi,Nx,deltax);
d_fourier_phi = 2i*pi/L*lin_n.*fourier_phi;
dphi = inverse_ft(d_fourier_phi,Nx,deltaf);

kin = sq(dphi)./2;
pot = (1/2*X.*X + 1/2*Beta*sq(phi) + V.*sin(k_scale.*X).^2).*sq(phi);
tot = kin + pot;
tot_L = tot(1:Nx-1);
tot_H = tot(2:Nx);
tot_avg = (tot_H+tot_L)/2;
y = sum(tot_avg)*deltax;
end

