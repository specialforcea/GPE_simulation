function y = integr2d(phi,Nx,Ny,deltax,deltay)
phi_L = phi(:,1:Ny-1);
phi_H = phi(:,2:Ny);
phi_avg = (phi_H+phi_L)/2;
y1 = sum(phi_avg,2)*deltay;
y = integr(y1,Nx,deltax);
end


