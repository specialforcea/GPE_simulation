l = [0.6 0.7 0.8 0.9];
lattice_pot = E_r/14.8*l(1);
V0 = lattice_pot/quanta;
Groundstate;
save('600Hz2.mat','phi_0');
lattice_pot = E_r/14.8*l(2);
V0 = lattice_pot/quanta;
Groundstate;
save('700Hz2.mat','phi_0');

lattice_pot = E_r/14.8*l(3);
V0 = lattice_pot/quanta;
Groundstate;
save('800Hz2.mat','phi_0');
lattice_pot = E_r/14.8*l(4);
V0 = lattice_pot/quanta;
Groundstate;
save('900Hz2.mat','phi_0');