phi1 = load('simulation_results/10242018evolve_with_soc/phi_less1_kR.mat');
phi1 = phi1.phi_less1_kR;
phi2 = load('simulation_results/10242018evolve_with_soc/phi_3kR.mat');
phi2 = phi2.phi_3kR;
phi3 = load('simulation_results/10242018evolve_with_soc/phi_4to5_kR_8193.mat');
phi3 = phi3.phi_4to5_kR_8193;
phi4 = load('simulation_results/10242018evolve_with_soc/phi_8_kR.mat');
phi4 = phi4.phi_8_kR;

speckle = load('speckle bench test data/T_simu/inten_10.mat');
speckle = speckle.inten;
speckle = 20*speckle./1e5;%average intensity about 5 in simulation units, times 20 to be 60% of initial chemical potential

mom_evo = zeros(5,20,2,200);

phi = phi3;
Raman_coupling = [1,5000, 10000, 20000, 40000];
Omega = Raman_coupling/Dip_freq;
final_phi = zeros(5,20,3,Nx);

for o=1:5

    for j=1:20
        phi = phi3;
        speckle_row = speckle(3*j,300:Nx+300);
        for i=1:200
            phi_1 = dynamic(phi,1e-3,1e-5,c0,c2,Nx,speckle_row,0,0,k_scale,f,deltax,deltaf,L,Omega(o),xmin,xmax,k_R,detuning);

            fp1 = sq(fourier_transform(phi_1(1,:),Nx,deltax));
            fp2 = sq(fourier_transform(phi_1(2,:),Nx,deltax));



            mom_evo(o,j,1,i) = sqrt(integr(f.^2.*fp1,Nx,deltaf));
            mom_evo(o,j,2,i) = sqrt(integr(f.^2.*fp2,Nx,deltaf));
            
            
            subplot(211)
            plot(f./k_spacing,fp1,f./k_spacing,fp2)
            subplot(212)
            plot(X,sq(phi_1))
                    

            phi = phi_1;

        end
        final_phi(o,j,:,:) = reshape(phi,[1,1,3,Nx]);
    end
end

save('simulation_results/10242018evolve_with_soc/mom_evo.mat','mom_evo');
save('simulation_results/10242018evolve_with_soc/final_phi.mat','final_phi');