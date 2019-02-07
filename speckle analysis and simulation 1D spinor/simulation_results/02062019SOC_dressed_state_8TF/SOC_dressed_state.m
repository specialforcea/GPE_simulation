%at detuning zero, ramp up omega in 25ms and hold for 25ms, get dressed
%state.

phi = phi_0;

ks = [1.0,1.3,1.5,1.7,2.0,2.2];
dts = (2+ks).^2 - ks.^2 - 4;
savepath = 'simulation_results/02062019SOC_dressed_state_8TF/phi_1_25ms_';

Omega = linspace(0,1*detuning,1601);

draw = 0;
spin1 = zeros(6,320);
spin2 = zeros(6,320);
width = zeros(6,320);
pos = zeros(6,320);


for j=1:8
    
    phi = phi_0;
    
for i=1:1600
    phi_1 = dynamic(phi,0.001,1e-5,c0,c2,Nx,0,1,0,k_scale,f,deltax,deltaf,L,Omega(i+1),xmin,xmax,k_R,dts(j)*detuning);
   % phi_1 = dynamic(phi,0.001,1e-5,c0,c2,Nx,0,1,0,k_scale,f,deltax,deltaf,L,4*detuning,xmin,xmax,k_R,detuning);

    phi = phi_1;
    if draw==1 && mod(i,60)==0
        figure(1)
        plot(X,sq(phi))
%         figure(2)
%         plot(f./k_R*2*pi,sq(fourier_transform(phi,Nx,deltax)));
         drawnow;
    end
    save(strcat(savepath,num2str(dts(j)),'d_1E.mat'),'phi_1');
end


% for i=1:320
%     phi_dress = dynamic(phi,0.005,1e-5,c0,c2,Nx,0,1,0,k_scale,f,deltax,deltaf,L,1*detuning,xmin,xmax,k_R,detuning);
%     phi = phi_dress;
%     if draw==1 && mod(i,10)==0
%         plot(X,sq(phi))
%         drawnow;
%     end
%     save('simulation_results/02062019SOC_dressed_state_8TF/phi_dress25ms_25ms_1d_1E.mat','phi_dress');
% end




%     phi = phip;
%     dt = linspace(detuning,0,st(j));


%     for i=1:st(j)
%         phi_dress = dynamic(phi,0.005,1e-5,c0,c2,Nx,0,1,0,k_scale,f,deltax,deltaf,L,1*detuning,xmin,xmax,k_R,dt(i));
%         phi = phi_dress;
%         if draw==1 && mod(i,10)==0
%             plot(X,sq(phi))
%             drawnow;
%         end
%         %save('simulation_results/01182019SOC_dressed_state/phi_dress25ms_25ms_0d_1E.mat','phi_dress');
%     end


phi = phi_1.*[exp(1i.*ks(j).*k_R.*X);exp(1i.*ks(j).*k_R.*X);exp(1i.*ks(j).*k_R.*X)];


    
    for i=1:320
        phi_dress = dynamic(phi,0.005,1e-5,c0,c2,Nx,0,0,0,k_scale,f,deltax,deltaf,L,1*detuning,xmin,xmax,k_R,0);
        phi = phi_dress;
        spin1(j,i) = integr(sq(phi(1,:)),Nx,deltax);
        spin2(j,i) = integr(sq(phi(2,:)),Nx,deltax);
        m1 = integr(X.*sq(phi(1,:)),Nx,deltax)/spin1(j,i);
        m2 = integr(X.*sq(phi(2,:)),Nx,deltax)/spin2(j,i);
        w1 = sqrt(integr((X-m1).^2.*sq(phi(1,:)),Nx,deltax)/spin1(j,i));
        w2 = sqrt(integr((X-m2).^2.*sq(phi(2,:)),Nx,deltax)/spin2(j,i));
        width(j,i) = w1*spin1(j,i) + w2*spin2(j,i);
        pos(j,i) = m1*spin1(j,i) + m2*spin2(j,i);
        if draw==1 && mod(i,30)==0
            plot(X,sq(phi(1,:)),X,sq(phi(2,:)))
            drawnow;
        end
   
    end
    
            save('simulation_results/02062019SOC_dressed_state_8TF/spin125ms_25ms.mat','spin1');
            save('simulation_results/02062019SOC_dressed_state_8TF/spin225ms_25ms.mat','spin2');
            save('simulation_results/02062019SOC_dressed_state_8TF/width25ms_25ms.mat','width');
            save('simulation_results/02062019SOC_dressed_state_8TF/pos25ms_25ms.mat','pos');
end


