%elolve for 48ms
scale_parameters;
Groundstate;

i=8;
filepath = strcat('speckle bench test data/numerical_speckle/13/inten_', num2str(i),'.mat');
speckle = load(filepath);
speckle = speckle.inten;

speckle = speckle/1e6;%average intensity about 1 in simulation units

speckle = speckle*5;%make it 68% of chemical potential.

X_shift = X;
X_shift(1:2800) = X_shift(1:2800) + L;

save_path = 'simulation_results/05162019kick_and_evolve_cluster/';
mom_evo = zeros(21,20,320);
prof_evo = zeros(21,20,320);
mean_m = zeros(21,20,320);
mean_p = zeros(21,20,320);
final_phi = zeros(21,20,Nx);
final_phi_200 = zeros(21,20,Nx);

kicks = (0.2:0.1:2.2);

for mk=1:21
    phi_mk = phi_0;
    phi_mk(1,:) = phi_mk(1,:).*exp(1i*k_R*kicks(mk)*X);
    rand_row = randi([1 Nx],1,20);
    for j=1:20

        speckle_row = speckle(rand_row(j),:);

        phi = phi_mk;

        for i=1:320
            phi_1 = dynamic(phi,0.005,1e-5,c0,c2,Nx,speckle_row,0,0,k_scale,f,deltax,deltaf,L,Omega,xmin,xmax);

            fp = sq(fourier_transform(phi_1(1,:),Nx,deltax));


            mean_mom = integr(f.*fp,Nx,deltaf);
            mean_prof = integr(X_shift.*sq(phi_1(1,:)),Nx,deltax);
            mean_m(mk,j,i) = mean_mom;
            mean_p(mk,j,i) = mean_prof;
            mom_evo(mk,j,i) = sqrt(integr((f-mean_mom).^2.*fp,Nx,deltaf));
            prof_evo(mk,j,i) = sqrt(integr((X_shift-mean_prof).^2.*sq(phi_1(1,:)),Nx,deltax));



            phi = phi_1;
            if i==200
                final_phi_200(mk,j,:) = reshape(phi(1,:),[1,1,Nx]);
                save(strcat(save_path,'final_phi_200.mat'),'final_phi_200');
            end

            
        end

        final_phi(mk,j,:) = reshape(phi(1,:),[1,1,Nx]);

        save(strcat(save_path,'mom_evo.mat'),'mom_evo')
        save(strcat(save_path,'prof_evo.mat'),'prof_evo')
        save(strcat(save_path,'final_phi.mat'),'final_phi')
        save(strcat(save_path,'mean_m.mat'),'mean_m')
        save(strcat(save_path,'mean_p.mat'),'mean_p')
    end
end
