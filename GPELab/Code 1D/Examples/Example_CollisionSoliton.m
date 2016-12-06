%%% This file is an example of how to use GPELab (FFT version)

%% Ground state of a Gross-Pitaevskii equation with quadratic potential and cubic nonlinearity in 1D


clear all;
%-----------------------------------------------------------
% Setting the data
%-----------------------------------------------------------

%% Setting the method and geometry
Computation = 'Dynamic';
Ncomponents = 2;
Type = 'Relaxation';
Deltat = 1e-2;
Stop_time = 25;
Stop_crit = [];
Method = Method_Var1d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);
xmin = -40;
xmax = 40;
Nx = 2^11+1;
Geometry1D = Geometry1D_Var1d(xmin,xmax,Nx);

%% Setting the physical problem
Delta = 1;
Beta = 1;
Physics1D = Physics1D_Var1d(Method,Delta,Beta); 
Physics1D = Dispersion_Var1d(Method, Physics1D);


alpha_1 = 0.25;
alpha_2 = -0.1965;
Coupled_NL{1,1} = @(Phi,X) -alpha_1*abs(Phi{1}).^2 - (alpha_1+2*alpha_2)*abs(Phi{2}).^2;
Coupled_NL{1,2} = @(Phi,X) 0;
Coupled_NL{2,1} = @(Phi,X) 0;
Coupled_NL{2,2} = @(Phi,X) -alpha_1*abs(Phi{2}).^2 - (alpha_1+2*alpha_2)*abs(Phi{1}).^2;
Physics1D = Nonlinearity_Var1d(Method, Physics1D, Coupled_NL); 

%% Setting the initial data
c_l = 1.2; 
n_l = 0.03;
b_l = sqrt(n_l+c_l^2/4);
X_l = -15;
X = Geometry1D.X ; 
Phi_0{1} = sqrt(2*Delta/abs(alpha_1))*b_l*sech(b_l*(X-X_l)).*exp(1i*c_l*X/2+ n_l);

c_r = -0.5; 
n_r = 0.1;
b_r = sqrt(n_r+c_r^2/4);
X_r = 0;
X = Geometry1D.X ; 
Phi_0{2} = sqrt(2*Delta/abs(alpha_1))*b_r*sech(b_r*(X-X_r)).*exp(1i*c_r*X/2+ n_r);

%% Setting informations and outputs
Solution_save = 1;
Outputs_iterations = 10;
Output_function{1} = @(Phi,X,FFTX) Geometry1D.dx*sum(X.*abs(Phi).^2);
Output_name{1} = 'Position of the soliton';
Outputs = OutputsINI_Var1d(Method,Outputs_iterations,Solution_save,Output_function,Output_name);
Printing = 1;
Evo = 100;
Draw = 1;
Print = Print_Var1d(Printing,Evo,Draw);

%-----------------------------------------------------------
% Launching simulation
%-----------------------------------------------------------

[Phi, Outputs] = GPELab1d(Phi_0,Method,Geometry1D,Physics1D,Outputs,[],Print);

Draw_Timesolution1d(Outputs,Method,Geometry1D,Figure_Var1d)

figure(3)
Time = [0:0.1:Stop_time]; 
plot(Time, Outputs.User_defined_local{1,1});
xlabel('Time');
ylabel('Position of Psi_1');
figure(4)
plot(Time, Outputs.User_defined_local{2,1});
xlabel('Time');
ylabel('Position of Psi_2');