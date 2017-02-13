function y = norm1d(phi, Nx, DeltaX)
phi_1 = conj(phi).*phi;
phi_L = phi_1(1:Nx-1);
phi_H = phi_1(2:Nx);
phi_avg = (phi_H+phi_L)/2;
y = sum(phi_avg)*DeltaX;
y = sqrt(y);
end