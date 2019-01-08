Rb_Mass = 1.42*10^-25;% mass of Rb87 in Kg
hbar = 1.05*10^-34;
k_b = 1.38*10^-23;
N = 2e5;%atom number
Dip_freq = 2*pi*10;%dipole trap frequency in Hz
a_0 = sqrt(hbar/Rb_Mass/Dip_freq); %characteristic length of wavefunction in dipole trap
raman_wavelength = 790*10^-9; % raman beam wave length
k = 2*pi/raman_wavelength; %raman beam wave vector
a_B = 5.29*10^-11; %Bohr radius
a_s0 = 101.8*a_B; %0 channel scattering length
a_s2 = 100.4*a_B; %2 channel scattering length
pre_c0 = 4*pi*hbar^2*(a_s0 + 2*a_s2)/(3*Rb_Mass);%real c0
pre_c2 = 4*pi*hbar^2*(a_s2 - a_s0)/(3*Rb_Mass);%real c2
c0 = pre_c0*N/(hbar*Dip_freq*a_0^3)/6;
c2 = pre_c2*N/(hbar*Dip_freq*a_0^3)/6;

quanta = hbar*Dip_freq;
k_scale = k*a_0;%dimensionless k in the probelm
Raman_beams_angle = 10;%deg
k_R = k_scale*sin(Raman_beams_angle/2/180*pi);

E_r = (hbar)^2*(k_R/a_0)^2/(2*Rb_Mass);%recoil energy
order = 3;%number of diffracted orders considered.
%Beta = 3000; dimensionless interaction.
TF_radius = (3*c0/2)^(1/3);
xmin = -4*TF_radius;
xmax = 4*TF_radius;

L = xmax-xmin;
Nx = 2^13+1; %number of grids.
Ns = 1280;

speckle_mag = 46;
deltas = 4.8e-6/speckle_mag/a_0;
smin = -(Ns-1)/2*deltas;
smax = (Ns-1)/2*deltas;

deltax = (xmax-xmin)/(Nx-1);
deltaf = 1/Nx/deltax;
deltasf = 1/Ns/deltax; 
X = linspace(xmin,xmax,Nx);% x space
S = linspace(smin,smax,Ns);
f = ((0:1:Nx-1)-(Nx-1)/2).*deltaf;% momentum space
sf = ((0:1:Ns-1)-(Ns-1)/2).*deltasf;
Raman_beams_angle = 10;%deg
k_spacing = k_R/2/pi;%2*k_scale in f space. 
lattice_pot = 30*E_r;
V = lattice_pot/quanta;
real_Omega = 20000;%Raman coupling in Hz
Omega = real_Omega/Dip_freq;%rescaled coupling strength
detuning = E_r/(hbar*Dip_freq);%recoil energy is defined with 2k_l
