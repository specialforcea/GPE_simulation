overlap1 = zeros(1,8);

for m = 1:4
    phi_m = phi_0;
    phi_m(1,:) = phi_m(1,:).*exp(1i*k_R*kick1(m)*X);
    
    phii = phi_m(1,:);
    ol = zeros(1,10);
    
    for j=1:10
        phif = reshape(final_phi(m,j,:),[1,Nx]);
        
        a = sqrt(integr(sq(phii),Nx,deltax));
        
        c = sqrt(integr(sq(phif),Nx,deltax));
        
        ki = fourier_transform(phii,Nx,deltax);
        
        kf = fourier_transform(phif,Nx,deltax);
        
        O = 1 - JSD(sq(ki),sq(kf),Nx,deltaf);
        
        ol(j) = a^2*c^2*O;
        
    end
    overlap1(m) = mean(ol);
    
end
        