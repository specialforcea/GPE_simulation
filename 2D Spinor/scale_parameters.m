Rb_Mass = 1.42*10^-25;% mass of Rb87 in Kg
hbar = 1.05*10^-34;
k_b = 1.38*10^-23;
N = 2e5;%atom number
Dip_freq = 2*pi*65;%x dipole trap frequency in Hz
gamma = 0.62;%y freq over x freq
a_0 = sqrt(hbar/Rb_Mass/Dip_freq); %characteristic length of wavefunction in dipole trap
raman_wavelength = 790*10^-9; % raman beam wave length
k = 2*pi/raman_wavelength; %raman beam wave vector
a_B = 5.29*10^-11; %Bohr radius
a_s0 = 101.8*a_B; %0 channel scattering length
a_s2 = 100.4*a_B; %2 channel scattering length
pre_c0 = 4*pi*hbar^2*(a_s0 + 2*a_s2)/(3*Rb_Mass);%real c0
pre_c2 = 4*pi*hbar^2*(a_s2 - a_s0)/(3*Rb_Mass);%real c2
c0 = pre_c0*N/(hbar*Dip_freq*a_0^3);
c2 = pre_c2*N/(hbar*Dip_freq*a_0^3);
E_r = (hbar)^2*(2*k)^2/(2*Rb_Mass);%recoil energy
quanta = hbar*Dip_freq;
k_scale = k*a_0;%dimensionless k in the probelm
order = 3;%number of diffracted orders considered.

TF_radius = (4*c0/pi)^(1/4);%Thomas Fermi radius
xmin = -TF_radius*6/5;
xmax = TF_radius*6/5;
ymin = -TF_radius/gamma*6/5;
ymax = TF_radius/gamma*6/5;
Lx = xmax-xmin;
Ly = ymax-ymin;
Nx = 2^8+1; %number of grids.
Ny = 2^8+1;
deltax = (xmax-xmin)/(Nx-1);
deltay = (ymax-ymin)/(Ny-1);
deltafx = 1/Nx/deltax;
deltafy = 1/Ny/deltay;
unity = ones(Nx,1);
X = unity * linspace(xmin,xmax,Nx);% x space
Y = linspace(ymin,ymax,Ny)' * unity';% x space

fx = ((0:1:Nx-1)-(Nx-1)/2).*deltafx;% X momentum space
fy = ((0:1:Ny-1)-(Ny-1)/2).*deltafy;% Y momentum space
k_spacing = k_scale*(xmax-xmin)/pi;%2*k_scale in f space. 
lattice_pot = E_r/4/4;
V = lattice_pot/quanta;
real_Omega = 20000;%Raman coupling in Hz
Omega = real_Omega/Dip_freq;%rescaled coupling strength
[paritx,parity] = fourier_parity(Nx,Ny);
dispersion = fourier_dispersion(Nx,Ny,Lx,Ly);
detuning = E_r/(hbar*Dip_freq);%recoil energy is defined with 2k_l
trap = 1;

