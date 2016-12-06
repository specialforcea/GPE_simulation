function AnyonsInter_y = AnyonsOperatorBy_2d(Phi,FFTX,FFTY)
epsilon = 1e-20;
AnyonsInter_y = 1i*ifft2(fft2(Phi).*FFTX./(FFTX.^2+FFTY.^2+epsilon));
AnyonsInter_y(isnan(AnyonsInter_y)) = 0;
end