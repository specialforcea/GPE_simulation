%at detuning zero, ramp up omega in 50ms and hold for 25ms, get dressed
%state.
scale_parameters;
Groundstate;
ks = (1.2:0.1:3.2);
OmegaR = 0.5;
band = (-4:0.001:4);
y1 = (band+1).^2;
y2 = (band-1).^2;
dt = abs(y2-y1);
up = min(y1,y2) + 1/2*(dt+sqrt(dt.^2+OmegaR^2));
down = min(y1,y2) + 1/2*(dt-sqrt(dt.^2+OmegaR^2));

k_minima = sqrt((4*ks).^2./(OmegaR^2 + (4*ks).^2));
dtn = 4.*(ks - k_minima);

% ks = [0.2,0.4,0.6,0.8];
% dts = (2+ks).^2 - ks.^2 - 4;
savepath = 'simulation_results/05282019SOC_dressed_state_8TF13/phi_1_50ms_';

Omega = linspace(0,OmegaR*detuning,3201);

draw = 0;
% spin1_2 = zeros(4,320);
% spin2_2 = zeros(4,320);
% width_2 = zeros(4,320);
% pos_2 = zeros(4,320);
for j=1:21
phi = phi_0;
for i=1:3200
    phi_1 = dynamic(phi,0.001,1e-5,c0,c2,Nx,0,1,0,k_scale,f,deltax,deltaf,L,Omega(i+1),xmin,xmax,k_R,dtn(j)*detuning);
   % phi_1 = dynamic(phi,0.001,1e-5,c0,c2,Nx,0,1,0,k_scale,f,deltax,deltaf,L,4*detuning,xmin,xmax,k_R,detuning);

    phi = phi_1;
    if draw==1 && mod(i,60)==0
        figure(1)
        plot(X,sq(phi))        
         drawnow;
    end
    
end


kick = ks(j) - k_minima(j);
phi_1 = phi_1.*[exp(1i.*kick.*k_R.*X);exp(1i.*kick.*k_R.*X);exp(1i.*kick.*k_R.*X)];
save(strcat(savepath,num2str(ks(j)),'k_R_',num2str(OmegaR),'E.mat'),'phi_1');
% spin1 = zeros(1,320);
% draw = 0;
%     
% for i=1:320
%     phi_dress = dynamic(phi,0.005,1e-5,c0,c2,Nx,0,0,0,k_scale,f,deltax,deltaf,L,OmegaR*detuning,xmin,xmax,k_R,0);
%     phi = phi_dress;
%     spin1(i) = integr(sq(phi(1,:)),Nx,deltax);
% 
%     if draw==1 && mod(i,30)==0
%         plot(1:1:320,spin1)
%         drawnow;
%     end
%     
% end
%     
% save(strcat(savepath,num2str(ks(j)),'k_R_',num2str(OmegaR),'E_spin1.mat'),'spin1');
% save(strcat(savepath,num2str(ks(j)),'k_R_',num2str(OmegaR),'E_fphi.mat'),'phi');
end


