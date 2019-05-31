
scale_parameters;
ks = (1.2:0.1:3.2);
OmegaR = 0.7;
k_minima = sqrt((4*ks).^2./(OmegaR^2 + (4*ks).^2));


i=8;
filepath = strcat('speckle bench test data/numerical_speckle/13/inten_', num2str(i),'.mat');
speckle = load(filepath);
speckle = speckle.inten;

speckle = speckle/1e6;%average intensity about 1 in simulation units

speckle = speckle*5;%make it 68% of chemical potential.


wavepath = 'simulation_results/05012019SOC_dressed_state_8TF13/phi_1_50ms_';
savepath = 'simulation_results/05062019kick_evolve_with_soc_cluster/07_4';


mom_evo = zeros(21,20,2,320);
prof_evo = zeros(21,20,2,320);
mean_m = zeros(21,20,2,320);
mean_p = zeros(21,20,2,320);
spin1 = zeros(21,20,320);
final_phi = zeros(21,20,2,Nx);




for mk=13:16

    load(strcat(wavepath,num2str(ks(mk)),'k_R_',num2str(OmegaR),'E.mat'));
    
	
	kick = ks(mk) - k_minima(mk);
	
    
    
    phi_mk = phi_1.*[exp(1i.*kick.*k_R.*X);exp(1i.*kick.*k_R.*X);exp(1i.*kick.*k_R.*X)];
    
        
        
    rand_row = randi([1 Nx],1,20);
    for j=1:20
		'alive'
        speckle_row = speckle(rand_row(j),:);

        phi = phi_mk;

        for i=1:320
            phi_1 = dynamic(phi,0.005,1e-5,c0,c2,Nx,speckle_row,0,0,k_scale,f,deltax,deltaf,L,OmegaR*detuning,xmin,xmax,k_R,0);
            
            spin1(mk,j,i) = integr(sq(phi_1(1,:)),Nx,deltax);
            fp1 = sq(fourier_transform(phi_1(1,:),Nx,deltax));
            fp2 = sq(fourier_transform(phi_1(2,:),Nx,deltax));

            mean_mom1 = integr(f.*fp1,Nx,deltaf)/spin1(mk,j,i);
            mean_mom2 = integr(f.*fp2,Nx,deltaf)/(1-spin1(mk,j,i));
            
            mean_prof1 = integr(X.*sq(phi_1(1,:)),Nx,deltax)/spin1(mk,j,i);
            mean_prof2 = integr(X.*sq(phi_1(2,:)),Nx,deltax)/(1-spin1(mk,j,i));

            mom_evo(mk,j,1,i) = sqrt(integr((f-mean_mom1).^2.*fp1,Nx,deltaf)/spin1(mk,j,i));
            mom_evo(mk,j,2,i) = sqrt(integr((f-mean_mom2).^2.*fp2,Nx,deltaf)/(1-spin1(mk,j,i)));
            
            prof_evo(mk,j,1,i) = sqrt(integr((X-mean_prof1).^2.*sq(phi_1(1,:)),Nx,deltax)/spin1(mk,j,i));
            prof_evo(mk,j,2,i) = sqrt(integr((X-mean_prof2).^2.*sq(phi_1(2,:)),Nx,deltax)/(1-spin1(mk,j,i)));
            
            mean_m(mk,j,1,i) = mean_mom1;
            mean_m(mk,j,2,i) = mean_mom2;
            
            mean_p(mk,j,1,i) = mean_prof1;
            mean_p(mk,j,2,i) = mean_prof2;

            phi = phi_1;
            
            %save(strcat(savepath,'/phi_',num2str(mk),'_',num2str(j),'_',num2str(i),'.mat'),'phi_1')

        end
        
        final_phi(mk,j,:,:) = reshape(phi(1:2,:),[1,1,2,Nx]);
        save(strcat(savepath,'/mom_evo.mat'),'mom_evo')
        save(strcat(savepath,'/prof_evo.mat'),'prof_evo')
        save(strcat(savepath,'/final_phi.mat'),'final_phi')
        save(strcat(savepath,'/mean_m.mat'),'mean_m')
        save(strcat(savepath,'/mean_p.mat'),'mean_p')
        save(strcat(savepath,'/spin1.mat'),'spin1')
    end
end
