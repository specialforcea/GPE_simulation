scale_parameters;
Groundstate;

i=6;
filepath = strcat('speckle bench test data/numerical_speckle/13/inten_', num2str(i),'.mat');
speckle = load(filepath);
speckle = speckle.inten;

speckle = speckle/1e6;

speckle = speckle*240; %600 Hz


savepath = 'simulation_results/11302020pulse_raman_and_speckle/detuning';


mom_evo = zeros(20,2,120);
prof_evo = zeros(20,2,120);
mean_m = zeros(20,2,120);
mean_p = zeros(20,2,120);
spin1 = zeros(20,30);
final_phi = zeros(20,20,2,Nx);






OmegaR = 1.0;
rand_row = randi([1 Nx],1,20);

for j=1:20
    j
    speckle_row = speckle(rand_row(j),:);

    phi = phi_0;

    for i=1:120
        phi_1 = dynamic(phi,0.001,1e-5,c0,c2,Nx,speckle_row,0,0,k_scale,f,deltax,deltaf,L,OmegaR*detuning,xmin,xmax,k_R,-4*detuning);

        spin1(j,i) = integr(sq(phi_1(1,:)),Nx,deltax);
        fp1 = sq(fourier_transform(phi_1(1,:),Nx,deltax));
        fp2 = sq(fourier_transform(phi_1(2,:),Nx,deltax));

        mean_mom1 = integr(f.*fp1,Nx,deltaf)/spin1(j,i);
        mean_mom2 = integr(f.*fp2,Nx,deltaf)/(1-spin1(j,i));

        mean_prof1 = integr(X.*sq(phi_1(1,:)),Nx,deltax)/spin1(j,i);
        mean_prof2 = integr(X.*sq(phi_1(2,:)),Nx,deltax)/(1-spin1(j,i));

        mom_evo(j,1,i) = sqrt(integr((f-mean_mom1).^2.*fp1,Nx,deltaf)/spin1(j,i));
        mom_evo(j,2,i) = sqrt(integr((f-mean_mom2).^2.*fp2,Nx,deltaf)/(1-spin1(j,i)));

        prof_evo(j,1,i) = sqrt(integr((X-mean_prof1).^2.*sq(phi_1(1,:)),Nx,deltax)/spin1(j,i));
        prof_evo(j,2,i) = sqrt(integr((X-mean_prof2).^2.*sq(phi_1(2,:)),Nx,deltax)/(1-spin1(j,i)));

        mean_m(j,1,i) = mean_mom1;
        mean_m(j,2,i) = mean_mom2;

        mean_p(j,1,i) = mean_prof1;
        mean_p(j,2,i) = mean_prof2;

        phi = phi_1;

        if mod(i,6)==0
            final_phi(int32(i/3),j,:,:) = reshape(phi(1:2,:),[1,1,2,Nx]);
        end
    end



    
    save(strcat(savepath,'/mom_evo.mat'),'mom_evo')
    save(strcat(savepath,'/prof_evo.mat'),'prof_evo')
    save(strcat(savepath,'/final_phi.mat'),'final_phi')
    save(strcat(savepath,'/mean_m.mat'),'mean_m')
    save(strcat(savepath,'/mean_p.mat'),'mean_p')
    save(strcat(savepath,'/spin1.mat'),'spin1')
end


