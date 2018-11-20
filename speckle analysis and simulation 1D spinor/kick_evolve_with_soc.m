i=4;
filepath = strcat('speckle bench test data/numerical_speckle/inten_', num2str(i),'.mat');
speckle = load(filepath);
speckle = speckle.inten;

speckle = 20*speckle/1e4*3;%average intensity about 5 in simulation units

mom_evo = zeros(6,10,2,300);
prof_evo = zeros(6,10,2,300);
mean_m = zeros(6,10,2,300);
mean_p = zeros(6,10,2,300);
final_phi = zeros(6,10,2,Nx);



for mk=5:6
    phi_mk = phi_dress;
    phi_mk = phi_mk.*repmat(exp(-1i*k_R*mk*X),3,1);
    rand_row = randi([1 Nx],1,10);
    for j=1:10

        speckle_row = speckle(rand_row(j),:);

        phi = phi_mk;

        for i=1:300
            phi_1 = dynamic(phi,1e-2,1e-5,c0,c2,Nx,speckle_row,0,0,k_scale,f,deltax,deltaf,L,2*detuning,xmin,xmax,k_R,detuning);
            
            fp1 = sq(fourier_transform(phi_1(1,:),Nx,deltax));
            fp2 = sq(fourier_transform(phi_1(2,:),Nx,deltax));

            mean_mom1 = integr(f.*fp1,Nx,deltaf);
            mean_mom2 = integr(f.*fp2,Nx,deltaf);
            
            mean_prof1 = integr(X.*sq(phi_1(1,:)),Nx,deltax);
            mean_prof2 = integr(X.*sq(phi_1(2,:)),Nx,deltax);

            mom_evo(mk,j,1,i) = sqrt(integr((f-mean_mom1).^2.*fp1,Nx,deltaf));
            mom_evo(mk,j,2,i) = sqrt(integr((f-mean_mom2).^2.*fp2,Nx,deltaf));
            
            prof_evo(mk,j,1,i) = sqrt(integr((X-mean_prof1).^2.*sq(phi_1(1,:)),Nx,deltax));
            prof_evo(mk,j,2,i) = sqrt(integr((X-mean_prof2).^2.*sq(phi_1(2,:)),Nx,deltax));
            
            mean_m(mk,j,1,i) = mean_mom1;
            mean_m(mk,j,2,i) = mean_mom2;
            
            mean_p(mk,j,1,i) = mean_prof1;
            mean_p(mk,j,2,i) = mean_prof2;


            
            

                    %plot(f(1,3500:4500)./k_spacing,fp(1,3500:4500))


            phi = phi_1;

        end
        final_phi(mk,j,:,:) = reshape(phi(1:2,:),[1,1,2,Nx]);
        save('simulation_results/11202018kick_evolve_with_soc_lap/mom_evo.mat','mom_evo')
        save('simulation_results/11202018kick_evolve_with_soc_lap/prof_evo.mat','prof_evo')
        save('simulation_results/11202018kick_evolve_with_soc_lap/final_phi.mat','final_phi')
        save('simulation_results/11202018kick_evolve_with_soc_lap/mean_m.mat','mean_m')
        save('simulation_results/11202018kick_evolve_with_soc_lap/mean_p.mat','mean_p')
    end
end
