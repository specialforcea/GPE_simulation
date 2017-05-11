function y = inverse_ft3(phi,paritx,parity,paritz,deltafx,deltafy,deltafz,Nx,Ny,Nz)


phi_1 = paritz.*parity.*paritx.*phi;

y = ifft2(phi_1).*deltafx*deltafy*deltafz*Ny*Nx*Nz;
y = y.*paritx.*parity.*paritz;
end