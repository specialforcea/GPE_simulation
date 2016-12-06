function Dipolar2d = RydBergInteraction2d(R0,Phi,FFTX,FFTY)
q = R0*sqrt(FFTX.^2+FFTY.^2);
K = R0*(2*pi^2/3)*(exp(-q/2)./q).*(exp(-q/2)-2*sin(pi/6-sqrt(3)/2*q));
K(isinf(K)) = R0*(2*pi^2/3);
Dipolar2d = real(ifft2(K.*fft2(abs(Phi).^2)));
end
