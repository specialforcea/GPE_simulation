Rb_Mass = 1.42*10^-25;% mass of Rb87 in Kg
hbar = 1.05*10^-34;
N = 1e5;%atom number
Dip_freq = 2*pi*70;%dipole trap frequency in Hz
a_0 = 1.3*10^-6; %characteristic length of wavefunction in dipole trap
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
%Beta = 3000; dimensionless interaction.
TF_radius = (15*c0/pi/4)^(1/5);

xmin = -TF_radius*6/5;
xmax = TF_radius*6/5;
ymin = -TF_radius*6/5;
ymax = TF_radius*6/5;
zmin = -TF_radius*6/5;
zmax = TF_radius*6/5;

Lx = xmax-xmin;
Ly = ymax-ymin;
Lz = zmax-zmin;

Nx = 2^8+1; %number of grids.
Ny = 2^8+1;
Nz = 2^8+1;

deltax = (xmax-xmin)/(Nx-1);
deltay = (ymax-ymin)/(Ny-1);
deltaz = (zmax-zmin)/(Nz-1);

deltafx = 1/(xmax-xmin);
deltafy = 1/(ymax-ymin);
deltafz = 1/(zmax-zmin);

unitx = ones(Ny,1);
Xa = unitx*linspace(xmin,xmax,Nx);% x space
Ya = linspace(ymin,ymax,Ny)'*ones(1,Nx);% x space
Za = linspace(zmin,zmax,Nz);% z space
X = Xa;
Y = Ya;
Z = zeros(Ny,Nx,Nz);
for i = 1:Nz
	X(:,:,i) = Xa;
	Y(:,;,i) = Ya;
	Z(:,:,i) = Za(i)*ones(Ny,Nx);
end





fx = (0:1:Nx-1)-(Nx-1)/2;% X momentum space
fy = (0:1:Ny-1)-(Ny-1)/2;% X momentum space
fz = (0:1:Nz-1)-(Nz-1)/2;% X momentum space

k_spacing = k_scale*(xmax-xmin)/pi;%2*k_scale in f space. 
lattice_pot = 4*E_r;
V = lattice_pot/quanta;
real_Omega = 20000;%Raman coupling in Hz
Omega = real_Omega/Dip_freq;%rescaled coupling strength
[paritx,parity,paritz] = fourier_parity(Nx,Ny,Nz);
dispersion = fourier_dispersion(Nx,Ny,Nz,Lx,Ly,Lz);
detuning = E_r/(hbar*Dip_freq);%recoil energy is defined with 2k_l