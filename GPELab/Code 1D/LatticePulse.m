Computation = 'Dynamic';
Ncomponents = 1;
Type = 'Splitting';
Deltat = 1e-5;
t_ramp = 0.2;
Stop_time = t_ramp;
Method = Method_Var1d(Computation,Ncomponents, Type, Deltat, Stop_time);
Method.Splitting = 'Lie';
xmin = -15;
xmax = 15;

Nx = 2^11+1;

Geometry1D = Geometry1D_Var1d(xmin,xmax, Nx);
Delta = 0.5;
Beta = 3000;



v = 200;
kappa = v/t_ramp;
k = 10;

Physics1D = TimePotential_Var1d(Method,Physics1D,@(t,x)(v*sin(k*x).^2));

Physics1D = Potential_Var1d(Method, Physics1D,@(x)(1/2)*(x.^2));

Physics1D = Nonlinearity_Var1d(Method, Physics1D,@(phi,x) abs(phi).^2);

Outputs_iterations = 100;
Outputs_save = 1;
Outputs = OutputsINI_Var1d(Method,Outputs_iterations,Outputs_save);
Printing = 1;
Evo = 100;
Draw = 1;
Print = Print_Var1d(Printing,Evo,Draw);
[Phi_2,Outputs]= GPELab1d(Phi_1,Method,Geometry1D,Physics1D,Outputs,[],Print);