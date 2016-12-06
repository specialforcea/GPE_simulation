Computation = 'Ground';
Ncomponents = 1;
Type = 'BESP';
Deltat = 1e-3;
Stop_time = [];
Stop_crit = {'MaxNorm',1e-8};
Method = Method_Var1d(Computation,Ncomponents, Type, Deltat, Stop_time ,Stop_crit);
xmin = -15;
xmax = 15;

Nx = 2^11+1;

Geometry1D = Geometry1D_Var1d(xmin, xmax, Nx); 
Delta = 0.5;
Beta = 3000; 

Optical_Potential = @(x) (1/2)*(x.^2);

Physics1D = Physics1D_Var1d(Method,Delta,Beta);
Physics1D = Potential_Var1d(Method, Physics1D,Optical_Potential);
Physics1D = Nonlinearity_Var1d(Method, Physics1D,@(phi,x) abs(phi).^2);
InitialData_choice = 2 ;
Phi_0 = InitialData_Var1d(Method, Geometry1D, Physics1D,InitialData_choice);
Outputs = OutputsINI_Var1d(Method);
Printing = 1;
Evo = 15;
Draw = 1;
Print = Print_Var1d(Printing,Evo,Draw);
[Phi_1,Outputs]= GPELab1d(Phi_0,Method,Geometry1D,Physics1D,Outputs,[],Print);