function AnyonsNL = AnyonsNonlinearity2d(Phi,FFTX,FFTY,AnyonsInter)
AnyonsQuadx = AnyonsInter*AnyonsOperatorBx_2d(abs(Phi).^2,FFTX,FFTY);
AnyonsQuady = AnyonsInter*AnyonsOperatorBy_2d(abs(Phi).^2,FFTX,FFTY);
AnyonsNL1 = AnyonsQuadx.^2 + AnyonsQuady.^2;
AnyonsNL2 = -2*(AnyonsInter*AnyonsOperatorBx_2d(AnyonsQuadx.*abs(Phi).^2,FFTX,FFTY)+AnyonsInter*AnyonsOperatorBy_2d(AnyonsQuady.*abs(Phi).^2,FFTX,FFTY));
SuperCourantx =  imag(Phi.*conj(ifft2(fft2(Phi).*1i*FFTX)));
SuperCouranty =  imag(Phi.*conj(ifft2(fft2(Phi).*1i*FFTY)));
AnyonsNL3 = -2*(AnyonsInter*AnyonsOperatorBx_2d(SuperCourantx,FFTX,FFTY)+AnyonsInter*AnyonsOperatorBy_2d(SuperCouranty,FFTX,FFTY));
AnyonsNL = AnyonsNL1 + AnyonsNL2 + AnyonsNL3;
end