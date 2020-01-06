create = 0;
analysis = 1;

% Nx = 2^12+1;
% deltax = (0.025)/(Nx-1);
% deltaf = 1/Nx/deltax;
% k_R = 2*pi/532e-9*12.5/sqrt(12.5^2+30^2);


if create==1
    for j=9:9
        N_x = 2^j + 1;
        X = (-2^(j-1):1:2^(j-1));
        X = repmat(X,2^j+1,1);
        Y = (-2^(j-1):1:2^(j-1))';
        Y = repmat(Y,1,2^j+1);
        beam = exp(-(X.^2 + Y.^2 )/2/71.7);

        rand_field = exp(-2*pi*1i*rand(2^j+1));
        iris = sqrt(X.^2 + Y.^2)<153.6;
        
        
        %field = field(250:3850,250:3850);
        field = beam.*rand_field.*iris;
        inten = sq(fftshift(fft2(field,Nx,Nx)));

        save(strcat('speckle bench test data/numerical_speckle/13/gauss_speckle_14mmbeam_15mmiris','.mat'), 'inten')


    end
end

if analysis==1
    N_x = 150;
    interval = (1:1:N_x);
    kr = interval*deltaf/k_R*2*pi;
    N_grid = size(kr,2);
    pow_r = zeros(10,N_grid);

    rr = interval*deltax*a_0;
    corr_r = zeros(30,N_grid);
    cors = zeros(1,30);

    for i=12:12
        filepath = strcat('speckle bench test data/numerical_speckle/13/gauss_speckle_14mmbeam_6mmiris','.mat');
        speckle = load(filepath);
        speckle = speckle.inten;
        speckle = double(speckle);
        speckle = speckle - mean(mean(speckle));
        fsp = fftshift(fft2(speckle));
        fsp = sq(fsp);

        corr = abs(fftshift(ifft2(fsp)));
        [maxr,maxc] = find(fsp==max(max(fsp)));
        %[maxr,maxc] = find(corr==max(max(corr)));


        X_cord = repmat((1:1:Nx),Nx,1);
        Y_cord = repmat((1:1:Nx)',1,Nx);

        X_cord = X_cord - mean(maxc);
        Y_cord = Y_cord - mean(maxr);

        r_cord = sqrt(X_cord.^2 + Y_cord.^2);
        for l=1:N_grid
            k = interval(l);
            [cordr,cordc] = find(r_cord<=k & r_cord>(k-1));
            num = size(cordr,1);
            po = 0;
            co = 0;
            for j=1:num
                po  = po + fsp(cordr(j),cordc(j));
                co = co + corr(cordr(j),cordc(j));
            end
            po = po/num;
            co = co/num;
            pow_r(i,l) = po;
            corr_r(i,l) = co;
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
        disp = 10;
        f = fit(kr(disp:end)',pow_r(i,disp:end)'./max(pow_r(i,disp:end)),'a*(acos(min(x/b,1)) - min(x/b,1)*sqrt(1-(min(x/b,1))^2))','StartPoint',[1/1.6,1]);
        coe = coeffvalues(f);
        plot(f,kr(disp:end),pow_r(i,disp:end)./max(pow_r(i,disp:end)))
        title('fit linear PSD to data')
        xlabel('k/k_R in radial direction')
        ylabel('normalized power density')
        text(3,0.3,strcat('k_m_a_x = ',num2str(coe(2))))


    %fit correlation length
%     figure(i)
%     corr_r(i,:) = corr_r(i,:)./max(corr_r(i,:));
%     f = fit(rr',corr_r(i,:)','a*exp(-x^2/2/b^2)+c','StartPoint',[1.0, 2e-6,0.05]);
%     coe = coeffvalues(f);
% 
%     plot(f,rr,corr_r(i,:))
%     title('correlation function averaged in radial direction')
%     xlabel('x/m')
%     ylabel('correlation')
%     text(0.4e-5,0.6,strcat('fit a*exp(-x^2/2/b^2), b = ',num2str(coe(2))))
%     cors(i) = coe(2);
    end
end
