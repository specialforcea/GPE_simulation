phi1 = load('simulation_results/10242018evolve_with_soc/phi_less1_kR.mat');
phi1 = phi1.phi_less1_kR;
phi2 = load('simulation_results/10242018evolve_with_soc/phi_3kR.mat');
phi2 = phi2.phi_3kR;
phi3 = load('simulation_results/10242018evolve_with_soc/phi_4to5_kR_8193.mat');
phi3 = phi3.phi_4to5_kR_8193;
phi4 = load('simulation_results/10242018evolve_with_soc/phi_8_kR.mat');
phi4 = phi4.phi_8_kR;

speckle1 = load('speckle bench test data/T_simu/inten_10.mat');
speckle1 = speckle1.inten;
speckle1 = 20*speckle1./1e5;%average intensity about 5 in simulation units, times 20 to be 60% of initial chemical potential

speckle2 = load('speckle bench test data/T_simu/inten_12.mat');
speckle2 = speckle2.inten;
speckle2 = 20*speckle2./1e5;%average intensity about 5 in simulation units, times 20 to be 60% of initial chemical potential


mom_evo_10 = zeros(3,10,2,300);
prof_evo_10 = zeros(3,10,2,300);
final_phi_10 = zeros(3,10,3,Nx);

mom_evo_12 = zeros(3,10,2,300);
prof_evo_12 = zeros(3,10,2,300);
final_phi_12 = zeros(3,10,3,Nx);



Omega = [0.001, 1,8]*detuning;


for o=1:3

    for j=1:10
        phi = phi3;
        speckle_row = speckle1(3*j,300:Nx+300);
        for i=1:300
            phi_1 = dynamic(phi,1e-2,1e-5,c0,c2,Nx,speckle_row,0,0,k_scale,f,deltax,deltaf,L,Omega(o),xmin,xmax,k_R,detuning);

            fp1 = sq(fourier_transform(phi_1(1,:),Nx,deltax));
            fp2 = sq(fourier_transform(phi_1(2,:),Nx,deltax));



            mom_evo_10(o,j,1,i) = sqrt(integr(f.^2.*fp1,Nx,deltaf));
            mom_evo_10(o,j,2,i) = sqrt(integr(f.^2.*fp2,Nx,deltaf));
            
            prof_evo_10(o,j,1,i) = sqrt(integr(X.^2.*sq(phi_1(1,:)),Nx,deltax));
            prof_evo_10(o,j,2,i) = sqrt(integr(X.^2.*sq(phi_1(2,:)),Nx,deltax));
            
            subplot(211)
            plot(f./k_spacing,fp1,f./k_spacing,fp2)
            subplot(212)
            plot(X,sq(phi_1))
                    

            phi = phi_1;

        end
        final_phi_10(o,j,:,:) = reshape(phi,[1,1,3,Nx]);
        save('simulation_results/10252018evolve_with_soc/mom_evo_10.mat','mom_evo_10');
        save('simulation_results/10252018evolve_with_soc/final_phi_10.mat','final_phi_10');
        save('simulation_results/10252018evolve_with_soc/prof_evo_10.mat','prof_evo_10');
    end
end

for o=1:3

    for j=1:10
        phi = phi3;
        speckle_row = speckle2(3*j,300:Nx+300);
        for i=1:300
            phi_1 = dynamic(phi,1e-2,1e-5,c0,c2,Nx,speckle_row,0,0,k_scale,f,deltax,deltaf,L,Omega(o),xmin,xmax,k_R,detuning);

            fp1 = sq(fourier_transform(phi_1(1,:),Nx,deltax));
            fp2 = sq(fourier_transform(phi_1(2,:),Nx,deltax));



            mom_evo_12(o,j,1,i) = sqrt(integr(f.^2.*fp1,Nx,deltaf));
            mom_evo_12(o,j,2,i) = sqrt(integr(f.^2.*fp2,Nx,deltaf));
            
            prof_evo_12(o,j,1,i) = sqrt(integr(X.^2.*sq(phi_1(1,:)),Nx,deltax));
            prof_evo_12(o,j,2,i) = sqrt(integr(X.^2.*sq(phi_1(2,:)),Nx,deltax));
            
            subplot(211)
            plot(f./k_spacing,fp1,f./k_spacing,fp2)
            subplot(212)
            plot(X,sq(phi_1))
                    

            phi = phi_1;

        end
        final_phi_12(o,j,:,:) = reshape(phi,[1,1,3,Nx]);
        save('simulation_results/10252018evolve_with_soc/mom_evo_12.mat','mom_evo_12');
        save('simulation_results/10252018evolve_with_soc/final_phi_12.mat','final_phi_12');
        save('simulation_results/10252018evolve_with_soc/prof_evo_12.mat','prof_evo_12');
    end
end
