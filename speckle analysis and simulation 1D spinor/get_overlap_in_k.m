
overlap1 = zeros(1,18);
wavepath = 'simulation_results/02152019SOC_dressed_state_8TF13/phi_1_50ms_';
OmegaR = 0.7;

for mk=1:18
    load(strcat(wavepath,num2str(dts(mk)),'d_',num2str(OmegaR),'E.mat'));
    phi_mk = phi_1.*[exp(1i.*ks(mk).*k_R.*X);exp(1i.*ks(mk).*k_R.*X);exp(1i.*ks(mk).*k_R.*X)];
    
    ol = zeros(1,10);
    for j=1:10
        phiiup = phi_mk(1,:);
        phiidown = phi_mk(2,:);
        
        phifup = reshape(final_phi(mk,j,1,:),[1,Nx]);
        phifdown = reshape(final_phi(mk,j,2,:),[1,Nx]);
        a = sqrt(integr(sq(phiiup),Nx,deltax));
        b = sqrt(integr(sq(phiidown),Nx,deltax));
        c = sqrt(integr(sq(phifup),Nx,deltax));
        d = sqrt(integr(sq(phifdown),Nx,deltax));
        
        kiup = fourier_transform(phiiup,Nx,deltax);
        kidown = fourier_transform(phiidown,Nx,deltax);
        kfup = fourier_transform(phifup,Nx,deltax);
        kfdown = fourier_transform(phifdown,Nx,deltax);
        
        Oup = 1-JSD(sq(kiup),sq(kfup),Nx,deltaf);
        Odown = 1-JSD(sq(kidown),sq(kfdown),Nx,deltaf);
        ol(j) = a^2*c^2*Oup + b^2*d^2*Odown + 2*a*b*c*d*sqrt(Oup*Odown);
    end
    overlap1(mk) = mean(ol);
end
        
        
        