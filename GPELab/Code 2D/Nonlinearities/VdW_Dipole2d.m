function Dipolar2d = VdW_Dipole2d(BetaDipole,R0,Psi,X,Y)
K = fft(1./(R0.^6+(X.^2+Y.^2).^3));
Dipolar2d = BetaDipole*ifft2(K.*fft2(abs(Psi).^2));
end
