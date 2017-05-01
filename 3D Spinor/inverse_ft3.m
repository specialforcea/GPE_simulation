function y = inverse_ft2(phi,paritx,parity,deltafx,deltafy,Nx,Ny)


phi_1 = paritz.*parity.*paritx.*phi;

y = ifft2(phi_1).*deltafx*deltafy*deltafz*Ny*Nx*Nz;
y = y.*paritx.*parity.*paritz;
end