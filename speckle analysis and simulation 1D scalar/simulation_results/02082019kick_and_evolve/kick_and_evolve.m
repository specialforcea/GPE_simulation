

i=8;
filepath = strcat('speckle bench test data/numerical_speckle/13/inten_', num2str(i),'.mat');
speckle = load(filepath);
speckle = speckle.inten;

speckle = speckle/1e6;%average intensity about 1 in simulation units

speckle = speckle*5;%make it 68% of chemical potential.

save_path = 'C:\Experiments\simulation\GPE_yuchen\GPE_simulation\speckle analysis and simulation 1D scalar\simulation_results';
mom_evo = zeros(8,10,320);
prof_evo = zeros(8,10,320);
mean_m = zeros(8,10,320);
mean_p = zeros(8,10,320);
final_phi = zeros(8,10,Nx);

kicks = [0.2,0.5,0.8,1.2,1.5,1.8,2.1,2.5];

for mk=1:8
    phi_mk = phi_0;
    phi_mk(1,:) = phi_mk(1,:).*exp(1i*k_R*kicks(mk)*X);
    rand_row = randi([1 Nx],1,10);
    for j=1:10

        speckle_row = speckle(rand_row(j),:);

        phi = phi_mk;

        for i=1:320
            phi_1 = dynamic(phi,0.005,1e-5,c0,c2,Nx,speckle_row,0,0,k_scale,f,deltax,deltaf,L,Omega,xmin,xmax);

            fp = sq(fourier_transform(phi_1(1,:),Nx,deltax));


            mean_mom = integr(f.*fp,Nx,deltaf);
            mean_prof = integr(X.*sq(phi_1(1,:)),Nx,deltax);
            mean_m(mk,j,i) = mean_mom;
            mean_p(mk,j,i) = mean_prof;
            mom_evo(mk,j,i) = sqrt(integr((f-mean_mom).^2.*fp,Nx,deltaf));
            prof_evo(mk,j,i) = sqrt(integr((X-mean_prof).^2.*sq(phi_1(1,:)),Nx,deltax));

                    %plot(f(1,3500:4500)./k_spacing,fp(1,3500:4500))


            phi = phi_1;

            save(strcat(save_path,'\02082019kick_and_evolve\phi_',num2str(mk),'_',num2str(j),'_',num2str(i),'.mat'),'phi_1')
        end

        final_phi(mk,j,:) = reshape(phi(1,:),[1,1,Nx]);

        save(strcat(save_path,'\02082019kick_and_evolve\mom_evo.mat'),'mom_evo')
        save(strcat(save_path,'\02082019kick_and_evolve\prof_evo.mat'),'prof_evo')
        save(strcat(save_path,'\02082019kick_and_evolve\final_phi.mat'),'final_phi')
        save(strcat(save_path,'\02082019kick_and_evolve\mean_m.mat'),'mean_m')
        save(strcat(save_path,'\02082019kick_and_evolve\mean_p.mat'),'mean_p')
    end
end
