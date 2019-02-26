ks = [0.2,0.4,0.6,0.8,1.0,1.3,1.5,1.7,2.0,2.2];
dts = (2+ks).^2 - ks.^2 - 4;
wavepath = 'simulation_results/02072019SOC_dressed_state_8TF13/phi_1_50ms_';

spin1 = zeros(10,200);
spin2 = zeros(10,200);
width = zeros(10,200);
pos = zeros(10,200);

draw = 0;

for j=1:10

load(strcat(wavepath,num2str(dts(j)),'d_1E.mat'));

phi = phi_1.*[exp(1i.*ks(j).*k_R.*X);exp(1i.*ks(j).*k_R.*X);exp(1i.*ks(j).*k_R.*X)];


    
    for i=1:200
        phi_dress = dynamic(phi,0.005,1e-5,0,0,Nx,0,0,0,k_scale,f,deltax,deltaf,L,1*detuning,xmin,xmax,k_R,0);
        phi = phi_dress;
        spin1(j,i) = integr(sq(phi(1,:)),Nx,deltax);
        spin2(j,i) = 1 - spin1(j,i);
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
    
    save('simulation_results/02252019kick_evolve_non_int/spin1_15ms_50ms.mat','spin1');
    save('simulation_results/02252019kick_evolve_non_int/spin2_15ms_50ms.mat','spin2');
    save('simulation_results/02252019kick_evolve_non_int/width_15ms_50ms.mat','width');
    save('simulation_results/02252019kick_evolve_non_int/pos_5ms_50ms.mat','pos');
end