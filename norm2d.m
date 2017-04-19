function y = norm2d(phi, Nx, Ny, deltax,deltay)
phi_1 = sq(phi);
phi_1 = sum(phi_1,3);
phi_Lx = phi_1(1:Nx-1,:);
phi_Hx = phi_1(2:Nx,:);
phi_avgx = (phi_Hx+phi_Lx)/2;
sumx = sum(phi_avgx,1)*deltax;

phi_Ly = sumx(1:Ny-1);
phi_Hy = sumx(2:Ny);
phi_avgy = (phi_Hy+phi_Ly)/2;
y = sum(phi_avgy)*deltay;

y = sqrt(y);
end