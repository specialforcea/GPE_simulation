

i=4;
filepath = strcat('speckle bench test data/numerical_speckle/inten_', num2str(i),'.mat');
speckle = load(filepath);
speckle = speckle.inten;

speckle = 20*speckle/1e4*3;%average intensity about 5 in simulation units



mom_evo = zeros(3,10,300);
prof_evo = zeros(3,10,300);
mean_m = zeros(3,10,300);
mean_p = zeros(3,10,300);
final_phi = zeros(3,10,Nx);

kick = [10,20,40];

for mk=1:3
    phi_mk = phi_0;
    phi_mk(1,:) = phi_mk(1,:).*exp(-1i*k_R*kick(mk)*X);
    rand_row = randi([1 Nx],1,10);
    for j=1:10

        speckle_row = speckle(rand_row(j),:);

        phi = phi_mk;

        for i=1:300
            phi_1 = dynamic(phi,1e-2,1e-5,c0,c2,Nx,speckle_row,0,0,k_scale,f,deltax,deltaf,L,Omega,xmin,xmax);

            fp = sq(fourier_transform(phi_1(1,:),Nx,deltax));


            mean_mom = integr(f.*fp,Nx,deltaf);
            mean_prof = integr(X.*sq(phi_1(1,:)),Nx,deltax);
            mean_m(mk,j,i) = mean_mom;
            mean_p(mk,j,i) = mean_prof;
            mom_evo(mk,j,i) = sqrt(integr((f-mean_mom).^2.*fp,Nx,deltaf));
            prof_evo(mk,j,i) = sqrt(integr((X-mean_prof).^2.*sq(phi_1(1,:)),Nx,deltax));

                    %plot(f(1,3500:4500)./k_spacing,fp(1,3500:4500))


            phi = phi_1;

        end
        final_phi(mk,j,:) = reshape(phi(1,:),[1,1,Nx]);
        save('simulation_results/11192018kick_and_evolve/mom_evo.mat','mom_evo')
        save('simulation_results/11192018kick_and_evolve/prof_evo.mat','prof_evo')
        save('simulation_results/11192018kick_and_evolve/final_phi.mat','final_phi')
        
        save('simulation_results/11192018kick_and_evolve/mean_m.mat','mean_m')
        save('simulation_results/11192018kick_and_evolve/mean_p.mat','mean_p')
    end
end
