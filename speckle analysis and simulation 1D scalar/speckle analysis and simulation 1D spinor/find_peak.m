function fp = find_peak(phi,order,k_spacing,Nx,deltaf)
fp = zeros(1,2*order +1);
%fp_index = zeros(1,2*order +1);
freq = (-order:1:order).*k_spacing;
freq_low = int16(freq - k_spacing/10) + (Nx-1)/2;
freq_high = int16(freq + k_spacing/10) + (Nx-1)/2;
for i = 1:2*order + 1
    %[fp(i) ,fp_index(i)] = max(phi(freq_low(i):freq_high(i)));
    fp(i) = integr(phi(freq_low(i):freq_high(i)),size(phi(freq_low(i):freq_high(i)),2),deltaf);
end

end

    



