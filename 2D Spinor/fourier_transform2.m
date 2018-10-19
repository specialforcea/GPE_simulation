function y = fourier_transform2(phi,deltax,deltay)
 


% phi_1 = paritx.*phi;
% phi_2 = parity.*phi_1;
% y = fft2(phi_2).*deltax.*deltay;
% y = y.* paritx.*parity;

y = fftshift(fft2(phi)).*deltax.*deltay;
end