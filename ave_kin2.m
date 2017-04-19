function y = ave_kin2(phi, Nx,deltax,deltaf,L)
lin_n =  (0:1:Nx-1) - (Nx-1)/2;
fourier_phi = fourier_transform(phi,Nx,deltax);
d_fourier_phi = (2*pi/L*lin_n).^2.*sq(fourier_phi)/2;
y = integr(d_fourier_phi,Nx,deltaf);
end