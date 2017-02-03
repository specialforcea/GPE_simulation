function y = integr(phi,Nx,DeltaX)
phi_L = phi(1:Nx-1);
phi_H = phi(2:Nx);
phi_avg = (phi_H+phi_L)/2;
y = sum(phi_avg)*DeltaX;

end