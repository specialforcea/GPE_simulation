create = 1;
analysis = 1;

if create==1
    %j = 6;
            N0 = 2000;
            N1 = 8001;
        
        field = exp(-2*pi*1i*rand(N0));%.*aperture;

deltas0 = 25e-3/N1/a_0;

smin0 = -(N1-1)/2*deltas0;
smax0 = (N1-1)/2*deltas0;


deltasf0 = 1/N1/deltas0;

        
        %field = field(250:3850,250:3850);
        inten = sq(fftshift(fft2(field,N1,N1)));

        %save(strcat('speckle bench test data/numerical_speckle/13/inten_', num2str(j),'.mat'), 'inten')
        %save(strcat('speckle bench test data/numerical_speckle/13/inten_Ns', num2str(j),'.mat'), 'inten')


    
end

if analysis==1
    N_x = 200;
    
kr = (1:1:N_x)*deltasf0/k_R*2*pi;
pow_r = zeros(3,N_x);

rr = (1:1:N_x)*deltas0*a_0;
corr_r = zeros(3,N_x);
cors = zeros(1,8);

%i = j;
    i = 1;
    %filepath = strcat('speckle bench test data/numerical_speckle/13/inten_Ns', num2str(i),'.mat');
    %speckle = load(filepath);
    %speckle = speckle.inten;
    speckle = inten;
    speckle = double(speckle);
    speckle = speckle - mean(mean(speckle));
    speckle = double(speckle);
    
    %fsp = fourier_transform(speckle',Ns,deltas)';
    fsp = fftshift(fft2(speckle));
    fp = sq(fsp);
    
    corr = abs(fftshift(ifft2(fp)));
    [maxr,maxc] = find(fp==max(max(fp)));
    
    
    
    X_cord = repmat((1:1:N1),N1,1);
    Y_cord = repmat((1:1:N1)',1,N1);
    
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
    f = fit(kr(7:end)',pow_r(i,7:end)'./max(pow_r(i,7:end)),'a*(acos(min(x/b,1)) - min(x/b,1)*sqrt(1-(min(x/b,1))^2))','StartPoint',[1/1.3,1]);
    coe = coeffvalues(f);
    plot(f,kr(7:end),pow_r(i,7:end)./max(pow_r(i,7:end)))
    title('fit linear PSD to data')
    xlabel('k/k_R in radial direction')
    ylabel('normalized psd')
    text(0.4,0.5,strcat('k_m_a_x = ',num2str(coe(2))))


%fit correlation length
% figure(i)
% corr_r(i,:) = corr_r(i,:)./max(corr_r(i,:));
% f = fit(rr',corr_r(i,:)','a*exp(-x^2/2/b^2)+c','StartPoint',[1.0, 2e-6,0.05]);
% coe = coeffvalues(f);
% 
% plot(f,rr,corr_r(i,:))
% title('correlation function averaged in radial direction')
% xlabel('x/m')
% ylabel('correlation')
% text(0.4e-5,0.6,strcat('fit a*exp(-x^2/2/b^2), b = ',num2str(coe(2))))
% cors(i) = coe(2);
end


