%elolve for 1ms
scale_parameters;
Groundstate;

i=6;
filepath = strcat('speckle bench test data/numerical_speckle/13/inten_', num2str(i),'.mat');
speckle = load(filepath);
speckle = speckle.inten;

speckle = speckle/1e6;%average intensity about 1 in simulation units

speckle = speckle*0;%make it 200Hz of chemical potential.

save_path = 'simulation_results/09032020speckle_pulse_exp_condition_cluster/1800/';
mom_evo = zeros(20,60);
prof_evo = zeros(20,60);
mean_m = zeros(20,60);
mean_p = zeros(20,60);
final_phi = zeros(20,20,Nx);


phi_mk = phi_0;
rand_row = randi([1 Nx],1,20);
for j=1:20
    j
    speckle_row = speckle(rand_row(j),:);
    size(speckle_row)
    phi = phi_mk;

    for i=1:60
        phi_1 = dynamic(phi,0.001,1e-5,c0,c2,Nx,speckle_row,0,0,k_scale,f,deltax,deltaf,L,Omega,xmin,xmax);

        fp = sq(fourier_transform(phi_1(1,:),Nx,deltax));


        mean_mom = integr(f.*fp,Nx,deltaf);
        mean_prof = integr(X.*sq(phi_1(1,:)),Nx,deltax);
        mean_m(j,i) = mean_mom;
        mean_p(j,i) = mean_prof;
        mom_evo(j,i) = sqrt(integr((f-mean_mom).^2.*fp,Nx,deltaf));
        prof_evo(j,i) = sqrt(integr((X-mean_prof).^2.*sq(phi_1(1,:)),Nx,deltax));

                %plot(f(1,3500:4500)./k_spacing,fp(1,3500:4500))


        phi = phi_1;
        if mod(i,3)==0
            final_phi(j,int32(i/3),:) = reshape(phi(1,:),[1,1,Nx]);
        end

        %save(strcat('simulation_results/03182019kick_and_evolve/phi_',num2str(mk),'_',num2str(j),'_',num2str(i),'.mat'),'phi_1')
    end

    %final_phi(mk,j,:) = reshape(phi(1,:),[1,1,Nx]);

    save(strcat(save_path,'mom_evo.mat'),'mom_evo')
    save(strcat(save_path,'prof_evo.mat'),'prof_evo')
    save(strcat(save_path,'final_phi.mat'),'final_phi')
    save(strcat(save_path,'mean_m.mat'),'mean_m')
    save(strcat(save_path,'mean_p.mat'),'mean_p')
end

