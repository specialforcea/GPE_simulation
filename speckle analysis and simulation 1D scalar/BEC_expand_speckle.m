
speckle1 = load('speckle bench test data/T_simu/inten_10.mat');
speckle1 = speckle1.inten;
speckle1 = speckle1./1e5;%average intensity about 5 in simulation units

speckle2 = load('speckle bench test data/T_simu/inten_12.mat');
speckle2 = speckle2.inten;
speckle2 = speckle2./1e5;%average intensity about 5 in simulation units


fac = [10,20,40];

mom_evo_10 = zeros(3,10,300);
prof_evo_10 = zeros(3,10,300);
final_phi_10 = zeros(3,10,Nx);


for o = 1:3
    speckle_ = speckle1*fac(o);
    for j=1:10
        speckle_row = speckle_(3*j,300:Nx+300);
        %speckle_row = speckle_row -mean(speckle_row);
        phi = phi_0;
        
        for i=1:300
            phi_1 = dynamic(phi,1e-2,1e-5,c0,c2,Nx,speckle_row,0,0,k_scale,f,deltax,deltaf,L,Omega,xmin,xmax);

            fp = sq(fourier_transform(phi_1(1,:),Nx,deltax));



            mom_evo_10(o,j,i) = sqrt(integr(f.^2.*fp,Nx,deltaf));
            prof_evo_10(o,j,i) = sqrt(integr(X.^2.*sq(phi_1(1,:)),Nx,deltax));

                    %plot(f(1,3500:4500)./k_spacing,fp(1,3500:4500))
                    

            phi = phi_1;

        end
        final_phi_10(o,j,:) = reshape(phi(1,:),[1,1,Nx]);
        save('simulation_results/10262018BEC_expand_speckle/mom_evo_10.mat','mom_evo_10')
        save('simulation_results/10262018BEC_expand_speckle/prof_evo_10.mat','prof_evo_10')
        save('simulation_results/10262018BEC_expand_speckle/final_phi_10.mat','final_phi_10')
    end
end



mom_evo_12 = zeros(3,10,300);
prof_evo_12 = zeros(3,10,300);
final_phi_12 = zeros(3,10,Nx);


for o = 1:3
    speckle_ = speckle2*fac(o);
    for j=1:10
        speckle_row = speckle_(3*j,300:Nx+300);
        %speckle_row = speckle_row -mean(speckle_row);
        phi = phi_0;
        
        for i=1:300
            phi_1 = dynamic(phi,1e-2,1e-5,c0,c2,Nx,speckle_row,0,0,k_scale,f,deltax,deltaf,L,Omega,xmin,xmax);

            fp = sq(fourier_transform(phi_1(1,:),Nx,deltax));



            mom_evo_12(o,j,i) = sqrt(integr(f.^2.*fp,Nx,deltaf));
            prof_evo_12(o,j,i) = sqrt(integr(X.^2.*sq(phi_1(1,:)),Nx,deltax));

                    %plot(f(1,3500:4500)./k_spacing,fp(1,3500:4500))
                    

            phi = phi_1;

        end
        final_phi_12(o,j,:) = reshape(phi(1,:),[1,1,Nx]);
        save('simulation_results/10262018BEC_expand_speckle/mom_evo_12.mat','mom_evo_12')
        save('simulation_results/10262018BEC_expand_speckle/prof_evo_12.mat','prof_evo_12')
        save('simulation_results/10262018BEC_expand_speckle/final_phi_12.mat','final_phi_12')
    end
end
