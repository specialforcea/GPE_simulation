%%% This file is an example of how to use GPELab (FFT version)

%% GROUND STATE COMPUTATION WITH A ROTATING TERM AND COUPLED NONLINEARITIES


%-----------------------------------------------------------
% Setting the data
%-----------------------------------------------------------

%% Setting the method and geometry
Computation = 'Ground';
Ncomponents = 3;
Type = 'BESP';
Deltat = 1e-2;
Stop_time = [];
Stop_crit = {'MaxNorm',1e-6};
Method = Method_Var2d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);
xmin = -10;
xmax = 10;
ymin = -10;
ymax = 10;
Nx = 2^7+1;
Ny = 2^7+1;
Geometry2D = Geometry2D_Var2d(xmin,xmax,ymin,ymax,Nx,Ny);

%% Setting the physical problem
Delta = 0.5;
Kappa = 0;
Beta = 1;
Physics2D = Physics2D_Var2d(Method, Delta, Beta);

RashbaDispersion{1,1} = @(FFTX,FFTY) Delta*(FFTX.^2+FFTY.^2);
RashbaDispersion{1,2} = @(FFTX,FFTY) Kappa*FFTX - 1i*Kappa*FFTY;
RashbaDispersion{1,3} = @(FFTX,FFTY) 0;
RashbaDispersion{2,1} = @(FFTX,FFTY) - Kappa*FFTX - 1i*Kappa*FFTY;
RashbaDispersion{2,2} = @(FFTX,FFTY) Delta*(FFTX.^2+FFTY.^2);
RashbaDispersion{2,3} = @(FFTX,FFTY) Kappa*FFTX - 1i*Kappa*FFTY;
RashbaDispersion{3,1} = @(FFTX,FFTY) 0;
RashbaDispersion{3,2} = @(FFTX,FFTY) - Kappa*FFTX - 1i*Kappa*FFTY;
RashbaDispersion{3,3} = @(FFTX,FFTY) Delta*(FFTX.^2+FFTY.^2);

Physics2D = Dispersion_Var2d(Method, Physics2D,RashbaDispersion);

mu = 0;

CoupledPotential{1,1} = @(X,Y) (1/2)*(X.^2+Y.^2);
CoupledPotential{1,2} = @(X,Y) mu*(X-1i*Y);
CoupledPotential{1,3} = @(X,Y) 0;
CoupledPotential{2,1} = @(X,Y) mu*(X+1i*Y);
CoupledPotential{2,2} = @(X,Y) (1/2)*(X.^2+Y.^2);
CoupledPotential{2,3} = @(X,Y) mu*(X-1i*Y);
CoupledPotential{3,1} = @(X,Y) 0;
CoupledPotential{3,2} = @(X,Y) mu*(X+1i*Y);
CoupledPotential{3,3} = @(X,Y) (1/2)*(X.^2+Y.^2);

Physics2D = Potential_Var2d(Method, Physics2D, CoupledPotential);

c_0 = 100;
c_2 = 2;

CoupledSpinNL{1,1} = @(Phi,X,Y) (c_0+c_2)*(abs(Phi{1}).^2+abs(Phi{2}).^2) + (c_0-c_2)*abs(Phi{3}).^2;
CoupledSpinNL{1,2} = @(Phi,X,Y) c_2*Phi{2}.*conj(Phi{3});
CoupledSpinNL{1,3} = @(Phi,X,Y) 0;
CoupledSpinNL{2,1} = @(Phi,X,Y) c_2*conj(Phi{2}).*Phi{3};
CoupledSpinNL{2,2} = @(Phi,X,Y) (c_0+c_2)*(abs(Phi{1}).^2+abs(Phi{3}).^2) + c_0*abs(Phi{2}).^2;
CoupledSpinNL{2,3} = @(Phi,X,Y) c_2*conj(Phi{2}).*Phi{1};
CoupledSpinNL{3,1} = @(Phi,X,Y) 0;
CoupledSpinNL{3,2} = @(Phi,X,Y) c_2*Phi{2}.*conj(Phi{1});
CoupledSpinNL{3,3} = @(Phi,X,Y) (c_0+c_2)*(abs(Phi{2}).^2+abs(Phi{3}).^2) + (c_0-c_2)*abs(Phi{1}).^2;

Physics2D = Nonlinearity_Var2d(Method, Physics2D,CoupledSpinNL);

%% Setting the initial data
InitialData_Choice = 1;
Phi_0 = InitialData_Var2d(Method, Geometry2D, Physics2D, InitialData_Choice);

%% Setting informations and outputs
Outputs = OutputsINI_Var2d(Method);
Printing = 1;
Evo = 100;
Draw = 1;
Print = Print_Var2d(Printing,Evo,Draw);

%-----------------------------------------------------------
% Launching computation
%-----------------------------------------------------------

[Phi_1, Outputs] = GPELab2d(Phi_0,Method,Geometry2D,Physics2D,Outputs,[],Print);


%% Setting the method and geometry
Computation = 'Dynamic';
Ncomponents = 3;
Type = 'Relaxation';
Deltat = 1e-3;
Stop_time = 1;
Stop_crit = [];
Method = Method_Var2d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);

%% Setting the physical problem
Delta = 0.5;
Beta = 1;
Physics2D = Physics2D_Var2d(Method, Delta, Beta);
Physics2D = Dispersion_Var2d(Method, Physics2D,RashbaDispersion);

mu = -1/sqrt(2);

CoupledPotential{1,1} = @(X,Y) (1/2)*(X.^2+Y.^2);
CoupledPotential{1,2} = @(X,Y) mu*(X-1i*Y);
CoupledPotential{1,3} = @(X,Y) 0;
CoupledPotential{2,1} = @(X,Y) mu*(X+1i*Y);
CoupledPotential{2,2} = @(X,Y) (1/2)*(X.^2+Y.^2);
CoupledPotential{2,3} = @(X,Y) mu*(X-1i*Y);
CoupledPotential{3,1} = @(X,Y) 0;
CoupledPotential{3,2} = @(X,Y) mu*(X+1i*Y);
CoupledPotential{3,3} = @(X,Y) (1/2)*(X.^2+Y.^2);

Physics2D = Potential_Var2d(Method, Physics2D, CoupledPotential);
Physics2D = Nonlinearity_Var2d(Method, Physics2D,CoupledSpinNL);

%% Setting informations and outputs
Outputs = OutputsINI_Var2d(Method);
Printing = 1;
Evo = 10;
Draw = 1;
Print = Print_Var2d(Printing,Evo,Draw);

%-----------------------------------------------------------
% Launching simulation
%-----------------------------------------------------------

[Phi, Outputs] = GPELab2d(Phi_1,Method,Geometry2D,Physics2D,Outputs,[],Print);