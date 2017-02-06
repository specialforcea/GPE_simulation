function y = inverse_ft2(phi,parity,paritx,deltafx,deltafy,Nx,Ny)


phi_1 = paritx.*phi;
phi_2 = parity.*phi_1;
y = ifft2(phi_2).*deltafx*deltafy*Ny*Nx;
y = y.*paritx.*parity;