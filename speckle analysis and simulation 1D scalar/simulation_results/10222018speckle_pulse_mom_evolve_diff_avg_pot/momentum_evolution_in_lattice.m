% evo_momentum = zeros(135,200);
% phi = phi_0;
% k_N = int16(k_spacing/deltaf); 
% for i=1:200
%     phi_1 = dynamic(phi,2e-5,1e-6,c0,c2,Nx,0,1,V,k_scale,f,deltax,deltaf,L,Omega,xmin,xmax);
%    
%     fp = sq(fourier_transform(phi_1(1,:),Nx,deltax));
%     
%     for j=1:21
%         evo_momentum(6*j+1,i) = integr(fp((Nx-1)/2+(-11+j)*k_N-50:(Nx-1)/2+(-11+j)*k_N+50),101,deltaf);
%     end
%     plot(f,fp)
%     drawnow
%     
%     phi = phi_1;
%     
% end
% 


evo_momentum = zeros(5,20,200);
mom_equil_dist = zeros(5,20,Nx);
phi = phi_0;


speckle = load('speckle bench test data/T_simu/inten_12.mat');
speckle = speckle.inten;
speckle = speckle./1e5;%average intensity about 5 in simulation units

fac = [1,5,10,20,40];
for o = 1:5
    speckle1 = speckle*fac(o);

    
    for j=1:20
        speckle_row = speckle1(3*j,300:Nx+300);
        %speckle_row = speckle_row -mean(speckle_row);
        phi = phi_0;
        fp_temp = zeros(1,Nx);

        for i=1:200
            phi_1 = dynamic(phi,1e-4,1e-6,c0,c2,Nx,speckle_row,1,0,k_scale,f,deltax,deltaf,L,Omega,xmin,xmax);

            fp = sq(fourier_transform(phi_1(1,:),Nx,deltax));
            if i>100
                fp_temp = fp_temp + fp;
            end


            evo_momentum(o,j,i) = sqrt(integr(f.^2.*fp,Nx,deltaf));

%             plot(f,fp)
%             drawnow

            phi = phi_1;

        end
        fp_temp = fp_temp./100;
        mom_equil_dist(o,j,:) = reshape(fp_temp,[1,1,Nx]); 
    end
end

save('simulation_results/10222018speckle_pulse_mom_evolve_diff_avg_pot/evo_momentum.mat','evo_momentum');
save('simulation_results/10222018speckle_pulse_mom_evolve_diff_avg_pot/mom_equil_dist.mat','mom_equil_dist');