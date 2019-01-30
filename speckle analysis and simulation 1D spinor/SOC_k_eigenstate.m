phi = phi_e;
grad = zeros(1,3201);


grad(1:1601) = linspace(0,5.0,1601);
grad(1602:end) = linspace(5.0,0,1600);
draw = 1;

for i=1:3200
    phi_1 = dynamic(phi,0.001,1e-5,c0,c2,Nx,-grad(i+1)*X,1,0,k_scale,f,deltax,deltaf,L,1*detuning,xmin,xmax,k_R,detuning);
   % phi_1 = dynamic(phi,0.001,1e-5,c0,c2,Nx,0,1,0,k_scale,f,deltax,deltaf,L,4*detuning,xmin,xmax,k_R,detuning);

    phi = phi_1;
    if draw==1 && mod(i,60)==0
        figure(1)
        plot(X,sq(phi))
%         figure(2)
%         plot(f./k_R*2*pi,sq(fourier_transform(phi,Nx,deltax)));
         drawnow;
    end
    save('simulation_results/01282019SOC_k_eigenstate/phi_1_50ms_318ms_0d_1E.mat','phi_1');
end


spin1 = zeros(1,320);
spin2 = zeros(1,320);
width = zeros(1,320);
pos = zeros(1,320);
for i=1:320
    phi_dress = dynamic(phi,0.005,1e-5,c0,c2,Nx,0,0,0,k_scale,f,deltax,deltaf,L,1*detuning,xmin,xmax,k_R,detuning);
    phi = phi_dress;
    spin1(i) = integr(sq(phi(1,:)),Nx,deltax);
    spin2(i) = integr(sq(phi(2,:)),Nx,deltax);
    m1 = integr(X.*sq(phi(1,:)),Nx,deltax);
    m2 = integr(X.*sq(phi(2,:)),Nx,deltax);
    w1 = sqrt(integr((X-m1).^2.*sq(phi(1,:)),Nx,deltax));
    w2 = sqrt(integr((X-m2).^2.*sq(phi(2,:)),Nx,deltax));
    width(i) = w1*spin1(i) + w2*spin2(i);
    pos(i) = m1*spin1(i) + m2*spin2(i);
    if draw==1 && mod(i,10)==0
        plot(X,sq(phi(1,:)),X,sq(phi(2,:)))
        drawnow;
    end
        save('simulation_results/01282019SOC_k_eigenstate/spin125ms_50ms_318ms.mat','spin1');
        save('simulation_results/01282019SOC_k_eigenstate/spin225ms_50ms_318ms.mat','spin2');
        save('simulation_results/01282019SOC_k_eigenstate/width25ms_50ms_318ms.mat','width');
        save('simulation_results/01282019SOC_k_eigenstate/pos25ms_50ms_318ms.mat','pos');

end
