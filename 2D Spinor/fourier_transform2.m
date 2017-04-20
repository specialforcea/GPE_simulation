function y = fourier_transform2(phi,paritx,parity,deltax,deltay)
 


phi_1 = paritx.*phi;
phi_2 = parity.*phi_1;
y = fft2(phi_2).*deltax.*deltay;
y = y.* paritx.*parity;
end