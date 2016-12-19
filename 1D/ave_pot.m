function y = ave_pot(phi,Nx,DeltaX,X,V,k_scale,Beta)
pot = (1/2*X.*X + Beta*sq(phi) + V.*sin(k_scale.*X).^2).*sq(phi);
tot = pot;
tot_L = tot(1:Nx-1);
tot_H = tot(2:Nx);
tot_avg = (tot_H+tot_L)/2;
y = sum(tot_avg)*DeltaX;
end