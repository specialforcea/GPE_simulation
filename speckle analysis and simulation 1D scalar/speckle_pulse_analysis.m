
mom = zeros(10,300);
mean_mom = zeros(10,300);
potential = zeros(10,300);
w_err = zeros(1,10);
sf = zeros(10,51,Nx);
ff = zeros(10,51,Nx);

for i=1:10
    file_path = strcat('simulation_results/09032020speckle_pulse_exp_condition_cluster_6_measure_poten_long/',num2str((i-1)*200),'/mom_evo.mat');
    strct = load(file_path);
    mom_w = strct.mom_evo;
    mom(i,:) = mean(mom_w);
    w_err(i) = std(mom_w(:,300));
    
    file_path = strcat('simulation_results/09032020speckle_pulse_exp_condition_cluster_6_measure_poten_long/',num2str((i-1)*200),'/mean_m.mat');
    strct = load(file_path);
    mom_e = strct.mean_m;
    mean_mom(i,:) = mean(mom_e);
    
    file_path = strcat('simulation_results/09032020speckle_pulse_exp_condition_cluster_6_measure_poten_long/',num2str((i-1)*200),'/poten.mat');
    strct = load(file_path);
    poten = strct.poten;
    potential(i,:) = mean(poten);
    
%     file_path = strcat('simulation_results/09032020speckle_pulse_exp_condition_cluster_5_0.05_freevo/',num2str((i-1)*200),'/final_phi.mat');
%     strct = load(file_path);
%     final_phi = strct.final_phi;
%     sf(i,:,:) = mean(sq(final_phi));
%     sp = 0;
%     for j=1:20
%         sp = sp + sq(fourier_transform(reshape(final_phi(j,:,:),[51,Nx]),Nx,deltax));
%     end
%     ff(i,:,:) = reshape(sp/20,[1,51,Nx]);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%t = [(1:1:60)*0.001/Dip_freq (1:1:250)*0.005/Dip_freq+0.06/Dip_freq] ;
t = (1:1:120)*0.001/Dip_freq;
sig0 = mom(1:9,1:120)/k_R*2*pi*3.3;
wid = sqrt(20^2+sig0.^2*18^2);
plot(t*1000,wid)
ylabel('width of momentum distribution/k_r')
xlabel('t (ms)')
legend('0Hz','200Hz','400Hz','600Hz','800Hz','1000Hz','1200Hz','1400Hz','1600Hz')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fff = reshape(ff(3,3:10,:),[8,Nx]);
% fff = mean(fff,1);

% fff = reshape(ff(3,51,:),[1,Nx]);
% xx = (0:0.01:4);
% kc = 0.8;
% y = 0.005*(acos(min(xx/kc,1)) - min(xx/kc,1).*sqrt(1-(min(xx/kc,1)).^2));
% subplot(326)
% plot(f/k_R*2*pi,fff,[-flip(xx) xx],[flip(y) y])
% ylim([0.0,0.02])
% xlim([-3,3])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% poten_ = potential(:,1)*10;
% 
% %mean_width = mean(mom(:,300:end),2);
% mean_width = mom(:,310);
% size(poten_)
% size(mean_width)
% %plot(poten_,(mean_width/a_0*2*pi).^2*hbar/2/Rb_Mass/pi,poten_,poten_)
% %plot(poten_(2:end),(mean_width(2:end)/a_0*2*pi).^2*hbar/2/Rb_Mass/pi-(mean_width(1)/a_0*2*pi).^2*hbar/2/Rb_Mass/pi,poten_,poten_)
% errorbar(poten_(2:end),(mean_width(2:end)/a_0*2*pi).^2*hbar/2/Rb_Mass/pi-(mean_width(1)/a_0*2*pi).^2*hbar/2/Rb_Mass/pi,(w_err(2:end)/a_0*2*pi).^2*hbar/2/Rb_Mass/pi)
% hold on 
% plot(poten_,poten_)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% i=6;
% filepath = strcat('speckle bench test data/numerical_speckle/13/inten_', num2str(i),'.mat');
% speckle = load(filepath);
% speckle = speckle.inten;
% 
% speckle = speckle/1e6;%average intensity about 1 in simulation units
% 
% speckle = speckle*80;%make it 200Hz of chemical potential.
% 
% phi = phi_0(1,:);
% 
% arr = zeros(1,10);
% 
% for i=1:10
%     arr(i) = integr(sq(phi).*speckle(i,:),Nx,deltax);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



