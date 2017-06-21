function y = chem_pot(phi,X,Y,Nx,Ny,Beta,k_scale,deltax,deltay,deltafx,deltafy,V,Lx,Ly,paritx,parity,gamma)
lin_nx =  (0:1:Nx-1) - (Nx-1)/2;
lin_ny =  (0:1:Ny-1) - (Ny-1)/2;
fourier_phi = fourier_transform2(phi,paritx(:,:,1),parity(:,:,1),deltax,deltay);
d_fourier_phix = (2i*pi/Lx)*(ones(Ny,1)*lin_nx).*fourier_phi;
d_fourier_phiy = (2i*pi/Ly)*(lin_ny'*ones(1,Nx)).*fourier_phi;
dphix = inverse_ft2(d_fourier_phix,paritx(:,:,1),parity(:,:,1),deltafx,deltafy,Nx,Ny);
dphiy = inverse_ft2(d_fourier_phiy,paritx(:,:,1),parity(:,:,1),deltafx,deltafy,Nx,Ny);

kin = sq(dphix)./2 + sq(dphiy)./2;
pot = (1/2*X.*X + 1/2*gamma^2*Y.*Y + 1/2*Beta*sq(phi) + V.*sin(k_scale.*X).^2).*sq(phi);
tot = kin + pot;

y = integr2d(tot,Nx,Ny,deltax,deltay);
end

