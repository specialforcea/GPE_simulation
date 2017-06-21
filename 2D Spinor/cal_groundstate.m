l = [0.6 0.9 1.1 1.3];


lattice_pot = E_r/14.8*l(3);
V0 = lattice_pot/quanta;
Groundstate;
save('1100Hz.mat','phi_0');
lattice_pot = E_r/14.8*l(4);
V0 = lattice_pot/quanta;
Groundstate;
save('1300Hz.mat','phi_0');