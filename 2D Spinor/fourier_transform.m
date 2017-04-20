function y = fourier_transform(phi,Nx,deltax)
 
lin_1 = linspace(0,Nx-1,Nx);
parity = (-1).^lin_1;
phi_1 =parity.*phi;
y = fft(phi_1).*deltax;
y = y.* parity;
end