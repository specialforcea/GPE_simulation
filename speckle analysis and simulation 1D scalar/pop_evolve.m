evo_time = 0.1;
Deltat_count = 1e-4;
Deltat = 1e-6;
evoN = int16(evo_time/Deltat_count);
population = zeros(3,evoN);
phi_evo = phi_0;
for i = 1:evoN
    phi_evo = dynamic(phi_evo,Deltat_count,Deltat,c0,c2,Nx,V,k_scale,f,deltax,deltaf,L,Omega);
    population(1,i) = integr(sq(phi_evo(1,:)),Nx,deltax);
    population(2,i) = integr(sq(phi_evo(2,:)),Nx,deltax);
    population(3,i) = integr(sq(phi_evo(3,:)),Nx,deltax);
end
popN = linspace(1e-6,1e-4,evoN);
plot(popN,population')
