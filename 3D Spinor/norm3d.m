function y = norm2d(phi, Nx, Ny, Nz,deltax,deltay,deltaz)
phi_1 = sq(phi);
phi_1 = sum(phi_1,4);
phi_Lx = phi_1(1:Ny-1,:,:);
phi_Hx = phi_1(2:Ny,:,:);
phi_avgx = (phi_Hx+phi_Lx)/2;
sumx = sum(phi_avgx,1)*deltax;

phi_Ly = sumx(1:Nx-1,:);
phi_Hy = sumx(2:Nx,:);
phi_avgy = (phi_Hy+phi_Ly)/2;
sumy = sum(phi_avgy,1)*deltay;

phi_Lz = sumy(1:Nz-1);
phi_Hz = sumy(2:Nz);
phi_avgz = (phi_Hz+phi_Lz)/2;
y = sum(phi_avgz)*deltaz;

y = sqrt(y);
end