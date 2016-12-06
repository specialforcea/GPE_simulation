function Dipolar2d = DipoleDipole2d(RhoDipole,phi,x,y,fftx,ffty)
    square_xi = fftx.^2+ffty.^2;
    K = 2-3*sqrt(pi)*sqrt(square_xi).*exp(square_xi).*erfc(sqrt(square_xi));
    K(isnan(K)) = 1-1./(2*square_xi(isnan(K)))+3./(4*square_xi(isnan(K)).^2);
    K(isinf(K)) = 1-1./(2*square_xi(isinf(K)))+3./(4*square_xi(isinf(K)).^2);
    Dipolar2d = RhoDipole*ifft2(K.*fft2(abs(phi).^2));
end