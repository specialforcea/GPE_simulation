function fp = find_peak2d(phi,k_spacing,Nx,Ny,deltafx,deltafy,deltax,deltay,paritx,parity)
fftphi = fourier_transform2(phi(1:Nx-1,1:Ny-1,:),paritx,parity,deltax,deltay);
sqphi = sq(fftphi);
fp = zeros(1,2);
%fp_index = zeros(1,2*order +1);
%freq = (-order:1:order).*k_spacing;
freq_low1 = int16(-k_spacing - k_spacing/10) + (Nx-1)/2;
freq_high1 = int16(-k_spacing + k_spacing/10) + (Nx-1)/2;

freq_low2 = int16(0.0 - k_spacing/10) + (Nx-1)/2;
freq_high2 = int16(0.0 + k_spacing/10) + (Nx-1)/2;

%[fp(i) ,fp_index(i)] = max(phi(freq_low(i):freq_high(i)));
%fp(i) = integr(phi(freq_low(i):freq_high(i)),size(phi(freq_low(i):freq_high(i)),2),deltaf);
fp(1) = integr2d(sqphi(freq_low2:freq_high2,freq_low2:freq_high2,1),size(sqphi(freq_low2:freq_high2,freq_low2:freq_high2,1),1),size(sqphi(freq_low2:freq_high2,freq_low2:freq_high2,2)),deltafx,deltafy);
fp(2) = integr2d(sqphi(freq_low2:freq_high2,freq_low1:freq_high1,2),size(sqphi(freq_low2:freq_high2,freq_low1:freq_high1,2),1),size(sqphi(freq_low2:freq_high2,freq_low1:freq_high1,2)),deltafx,deltafy);

end

    



