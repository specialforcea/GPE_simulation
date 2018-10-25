
speckle = load('speckle bench test data/T_simu/inten_10.mat');
speckle = speckle.inten;
speckle = speckle./1e5;%average intensity about 5 in simulation units



fac = [1,5,10,20,40];

mom_evo = zeros(5,20,300);
prof_evo = zeros(5,20,300);
final_phi = zeros(5,20,Nx);

phi = phi_0;
for o = 1:5
    speckle1 = speckle*fac(o);
    for j=1:20
        speckle_row = speckle1(3*j,300:Nx+300);
        %speckle_row = speckle_row -mean(speckle_row);
        phi = phi_0;
        
        for i=1:300
            phi_1 = dynamic(phi,1e-2,1e-5,c0,c2,Nx,speckle_row,0,0,k_scale,f,deltax,deltaf,L,Omega,xmin,xmax);

            fp = sq(fourier_transform(phi_1(1,:),Nx,deltax));



            mom_evo(1,j,i) = sqrt(integr(f.^2.*fp,Nx,deltaf));
            prof_evo(1,j,i) = sqrt(integr(X.^2.*sq(phi_1(1,:)),Nx,deltax));

                    plot(f(1,3500:4500)./k_spacing,fp(1,3500:4500))
                    

            phi = phi_1;

        end
        final_phi(o,j,:) = reshape(phi(1,:),[1,1,Nx]);
        save('simulation_results/10252018BEC_expand_speckle/mom_evo.mat','mom_evo')
        save('simulation_results/10252018BEC_expand_speckle/prof_evo.mat','prof_evo')
        save('simulation_results/10252018BEC_expand_speckle/final_phi.mat','final_phi')
    end
end


