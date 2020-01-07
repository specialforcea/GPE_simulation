
draw=0;

spin1 = zeros(6,320);
spin2 = zeros(6,320);
width = zeros(6,320);
kick = linspace(0,2,6);

for j=1:1
    phi = phi_1;
    

    for i=1:320
        phi_dress = dynamic(phi,0.005,1e-5,c0,c2,Nx,0,0,0,k_scale,f,deltax,deltaf,L,1*detuning,xmin,xmax,k_R,detuning);
        phi = phi_dress;
        spin1(j,i) = integr(sq(phi(1,:)),Nx,deltax);
        spin2(j,i) = integr(sq(phi(2,:)),Nx,deltax);
        m1 = integr(X.*sq(phi(1,:)),Nx,deltax);
        m2 = integr(X.*sq(phi(2,:)),Nx,deltax);
        w1 = sqrt(integr((X-m1).^2.*sq(phi(1,:)),Nx,deltax));
        w2 = sqrt(integr((X-m2).^2.*sq(phi(2,:)),Nx,deltax));
        width(j,i) = w1*spin1(j,i) + w2*spin2(j,i);
        if draw==1 && mod(i,20)==0
            plot(X,sq(phi))
            drawnow;
        end
%         save('simulation_results/01222019SOC_dressed_state_8TF/spin125ms_318ms_free.mat','spin1');
%         save('simulation_results/01222019SOC_dressed_state_8TF/spin225ms_318ms_free.mat','spin2');
%         save('simulation_results/01222019SOC_dressed_state_8TF/width25ms_318ms_free.mat','width');
    end
end