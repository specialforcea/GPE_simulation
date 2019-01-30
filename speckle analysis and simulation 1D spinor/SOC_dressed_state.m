%at detuning zero, ramp up omega in 25ms and hold for 25ms, get dressed
%state.

% phi = phi_0;
% grad = zeros(1,1601);
% 
% Omega = linspace(0,1*detuning,1601);
% % grad(1:8001) = linspace(0,0.1,8001);
% % grad(8002:end) = linspace(0.1,0,8000);
% draw = 1;
% 
% for i=1:1600
%     phi_1 = dynamic(phi,0.001,1e-5,c0,c2,Nx,0,1,0,k_scale,f,deltax,deltaf,L,Omega(i+1),xmin,xmax,k_R,detuning);
%    % phi_1 = dynamic(phi,0.001,1e-5,c0,c2,Nx,0,1,0,k_scale,f,deltax,deltaf,L,4*detuning,xmin,xmax,k_R,detuning);
% 
%     phi = phi_1;
%     if draw==1 && mod(i,60)==0
%         figure(1)
%         plot(X,sq(phi))
% %         figure(2)
% %         plot(f./k_R*2*pi,sq(fourier_transform(phi,Nx,deltax)));
%          drawnow;
%     end
%     save('simulation_results/01282019SOC_dressed_state_30deg/phi_1_25ms_0d_1E.mat','phi_1');
% end


% for i=1:3200
%     phi_dress = dynamic(phi,0.005,1e-5,c0,c2,Nx,0,1,0,k_scale,f,deltax,deltaf,L,1*detuning,xmin,xmax,k_R,detuning);
%     phi = phi_dress;
%     if draw==1 && mod(i,10)==0
%         plot(X,sq(phi))
%         drawnow;
%     end
%     %save('simulation_results/01182019SOC_dressed_state/phi_dress25ms_25ms_0d_1E.mat','phi_dress');
% end
% 
spin1 = zeros(1,3200);
spin2 = zeros(1,3200);
width = zeros(1,3200);
pos = zeros(1,3200);
for i=1:3200
    phi_dress = dynamic(phi,0.005,1e-5,c0,c2,Nx,0,1,0,k_scale,f,deltax,deltaf,L,1*detuning,xmin,xmax,k_R,detuning);
    phi = phi_dress;
    spin1(i) = integr(sq(phi(1,:)),Nx,deltax);
    spin2(i) = integr(sq(phi(2,:)),Nx,deltax);
    m1 = integr(X.*sq(phi(1,:)),Nx,deltax);
    m2 = integr(X.*sq(phi(2,:)),Nx,deltax);
    w1 = sqrt(integr((X-m1).^2.*sq(phi(1,:)),Nx,deltax));
    w2 = sqrt(integr((X-m2).^2.*sq(phi(2,:)),Nx,deltax));
    width(i) = w1*spin1(i) + w2*spin2(i);
    pos(i) = m1*spin1(i) + m2*spin2(i);
    if draw==1 && mod(i,30)==0
        plot(X,sq(phi(1,:)),X,sq(phi(2,:)))
        drawnow;
    end
%             save('simulation_results/01242019SOC_dressed_state_8TF/spin125ms_250ms_grad01_free.mat','spin1');
%         save('simulation_results/01242019SOC_dressed_state_8TF/spin225ms_250ms_grad01_free.mat','spin2');
%         save('simulation_results/01242019SOC_dressed_state_8TF/width25ms_250ms_grad01_free.mat','width');
end


