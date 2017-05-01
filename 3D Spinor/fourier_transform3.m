function y = fourier_transform3(phi,paritx,parity,paritz,deltax,deltay,deltaz)
 


phi_1 = paritz.*parity.*paritx.*phi;

y = fftn(phi_1).*deltax.*deltay.*deltaz;
y = y.* paritx.*parity.*paritz;
end