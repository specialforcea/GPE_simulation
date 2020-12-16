%at detuning zero, ramp up omega in 50ms and hold for 16ms, get dressed
%state.
scale_parameters;
Groundstate;

savepath = 'simulation_results/09172020SOC_dressed_state_minima_100ms_cluster/phi_1_100ms_';

Omegas = [2.5,3.0,3.5,4.0];

for j=1:4
j
phi = phi_0;
OmegaR = Omegas(j);
band = (-4:0.001:4);
Deltaq = 4.*band;
y1 = (band+1).^2;
y2 = (band-1).^2;
up = Deltaq./2 + 1/2*sqrt(Deltaq.^2+OmegaR^2) + (band-1).^2;
down = Deltaq./2 - 1/2*sqrt(Deltaq.^2+OmegaR^2) + (band-1).^2;

Omega = linspace(0,OmegaR*detuning,6401);

for i=1:6400
    phi_1 = dynamic(phi,0.001,1e-5,c0,c2,Nx,0,1,0,k_scale,f,deltax,deltaf,L,Omega(i+1),xmin,xmax,k_R,0);
   
    phi = phi_1;
    
end

save(strcat(savepath,num2str(OmegaR),'E.mat'),'phi_1');

spin1 = zeros(1,200);
    
for i=1:200
    phi_dress = dynamic(phi,0.005,1e-5,c0,c2,Nx,0,0,0,k_scale,f,deltax,deltaf,L,OmegaR*detuning,xmin,xmax,k_R,0);
    phi = phi_dress;
    spin1(i) = integr(sq(phi(1,:)),Nx,deltax);

end
    
save(strcat(savepath,num2str(OmegaR),'E_spin1.mat'),'spin1');
end


