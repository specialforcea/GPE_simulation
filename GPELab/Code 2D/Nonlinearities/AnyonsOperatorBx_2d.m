function AnyonsInter_x = AnyonsOperatorBx_2d(Phi,FFTX,FFTY)
epsilon = 1e-20;
AnyonsInter_x = -1i*ifft2(fft2(Phi).*FFTY./(FFTX.^2+FFTY.^2+epsilon));
AnyonsInter_x(isnan(AnyonsInter_x)) = 0;
end