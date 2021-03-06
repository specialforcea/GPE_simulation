%elolve for 1ms
scale_parameters;
Groundstate;

i=6;
filepath = strcat('speckle bench test data/numerical_speckle/13/inten_', num2str(i),'.mat');
speckle = load(filepath);
speckle = speckle.inten;

speckle = speckle/1e6;


mom_evo = zeros(20,253);
prof_evo = zeros(20,253);
mean_m = zeros(20,253);
mean_p = zeros(20,253);
poten = zeros(20,3);
final_phi = zeros(20,51,Nx);

for k=1:10
    speckle_ = speckle*(k-1)*80;
    save_path = strcat('simulation_results/09032020speckle_pulse_exp_condition_cluster_6_0.05_freevo/',num2str((k-1)*200),'/');
    
    phi_mk = phi_0;
    rand_row = randi([1 Nx],1,20);
    k
    for j=1:20
    
        j
        speckle_row = speckle_(rand_row(j),:);

        phi = phi_mk;

        for i=1:3
            
            poten(j,i) = integr(sq(phi(1,:)).*speckle_row,Nx,deltax);
            
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
        
        
        for i=4:253
           
            
            phi_1 = dynamic(phi,0.005,1e-5,c0,c2,Nx,0,0,0,k_scale,f,deltax,deltaf,L,Omega,xmin,xmax);

            fp = sq(fourier_transform(phi_1(1,:),Nx,deltax));


            mean_mom = integr(f.*fp,Nx,deltaf);
            mean_prof = integr(X.*sq(phi_1(1,:)),Nx,deltax);
            mean_m(j,i) = mean_mom;
            mean_p(j,i) = mean_prof;
            mom_evo(j,i) = sqrt(integr((f-mean_mom).^2.*fp,Nx,deltaf));
            prof_evo(j,i) = sqrt(integr((X-mean_prof).^2.*sq(phi_1(1,:)),Nx,deltax));

                    %plot(f(1,3500:4500)./k_spacing,fp(1,3500:4500))


            phi = phi_1;
            if mod(i-3,5)==0
                final_phi(j,int32((i+2)/5),:) = reshape(phi(1,:),[1,1,Nx]);
            end

            %save(strcat('simulation_results/03182019kick_and_evolve/phi_',num2str(mk),'_',num2str(j),'_',num2str(i),'.mat'),'phi_1')
        end

        %final_phi(mk,j,:) = reshape(phi(1,:),[1,1,Nx]);

        save(strcat(save_path,'mom_evo.mat'),'mom_evo')
        save(strcat(save_path,'prof_evo.mat'),'prof_evo')
        save(strcat(save_path,'final_phi.mat'),'final_phi')
        save(strcat(save_path,'mean_m.mat'),'mean_m')
        save(strcat(save_path,'mean_p.mat'),'mean_p')
        save(strcat(save_path,'poten.mat'),'poten')
    end
end
