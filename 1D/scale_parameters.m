Rb_Mass = 1.42*10^-25;% mass of Rb87 in Kg
hbar = 1.05*10^-34;
Dip_freq = 2*pi*70;%dipole trap frequency in Hz
a_0 = 1.3*10^-6; %characteristic length of wavefunction in dipole trap
raman_wavelength = 790*10^-9; % raman beam wave length
k = 2*pi/raman_wavelength; %raman beam wave vector
a_B = 5.29*10^-11; %Bohr radius
E_r = (hbar)^2*(2*k)^2/(2*Rb_Mass);%recoil energy
quanta = hbar*Dip_freq;
k_scale = k*a_0;%dimensionless k in the probelm
order = 2;%number of diffracted orders considered.
Beta = 3000; %dimensionless interaction.
TF_radius = (3*Beta/2)^(1/3);
xmin = -TF_radius;
xmax = TF_radius;
Nx = 2^11+1; %number of grids.
DeltaX = (xmax-xmin)/(Nx-1);
X = linspace(xmin,xmax,Nx);% x space
f = (0:1:Nx-1)-(Nx-1)/2;% momentum space
k_spacing = k_scale*(xmax-xmin)/pi;%2*k_scale in f space. 
