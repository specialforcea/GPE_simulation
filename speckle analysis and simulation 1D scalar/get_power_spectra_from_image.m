li=[20,25,30];

N_x = 300;
kr = (1:1:N_x)*deltasf/k_R*4*pi;
pow_r = zeros(3,N_x);

rr = (1:1:N_x)*4.8e-6/46;
corr_r = zeros(3,N_x);


for i=1:3
    speckle = imread(strcat('speckle bench test data/12/',num2str(li(i)),'mm.png'));
    speckle = double(speckle);
    speckle = speckle - repmat(mean(speckle,2),1,Ns);
    
    %fsp = fourier_transform(speckle',Ns,deltas)';
    fsp = fftshift(fft2(speckle));
    fp = sq(fsp);
    
    corr = abs(fftshift(ifft2(fp)));
    [maxr,maxc] = find(fp==max(max(fp)));
    
    
    
    X_cord = repmat((1:1:1280),1024,1);
    Y_cord = repmat((1:1:1024)',1,1280);
    
    X_cord = X_cord - mean(maxc);
    Y_cord = Y_cord - mean(maxr);
    
    r_cord = sqrt(X_cord.^2 + Y_cord.^2);
    for k=1:N_x
        [cordr,cordc] = find(r_cord<=k & r_cord>(k-1));
        num = size(cordr,1);
        po = 0;
        co = 0;
        for j=1:num
            po  = po + fp(cordr(j),cordc(j));
            co = co + corr(cordr(j),cordc(j));
        end
        po = po/num;
        co = co/num;
        pow_r(i,k) = po;
        corr_r(i,k) = co;
    end
        
%plot PSD in linear and log scale%%        
%     figure(i)
%     subplot(121)
%     plot(r(3:end),pow_r(i,3:end))
%     title('PDS')
%     xlabel('k/k_R in radial direction')
%     ylabel('power density')
%     subplot(122)
%     plot(r,log(pow_r(i,:)))
%     title('log PDS')
%     xlabel('k/k_R in radial direction')
%     ylabel('log power density')
  

%fit linear PSD to data
    figure(i)
    f = fit(kr(7:end)',pow_r(i,7:end)','a*max(1-x/b,0)','StartPoint',[max(pow_r(i,7:end)),4]);
    coe = coeffvalues(f);
    plot(f,kr(7:end),pow_r(i,7:end))
    title('fit linear PSD to data')
    xlabel('k/k_R in radial direction')
    ylabel('power density')
    text(5,2e11,strcat('k_m_a_x = ',num2str(coe(2))))


%fit correlation length
% figure(i)
% corr_r(i,:) = corr_r(i,:)./max(corr_r(i,:));
% f = fit(rr',corr_r(i,:)','a*exp(-x^2/2/b^2)','StartPoint',[1.0, 2e-6]);
% coe = coeffvalues(f);
% 
% plot(f,rr,corr_r(i,:))
% title('correlation function averaged in radial direction')
% xlabel('x/m')
% ylabel('correlation')
% text(0.4e-5,0.6,strcat('fit a*exp(-x^2/2/b^2), b = ',num2str(coe(2))))

end
    
    
    