
Computation = 'Dynamic';
Ncomponents = 1;
Type = 'Splitting';
Deltat = 1e-3;
Stop_time = 1;
Stop_crit = [];
Method = Method_Var2d(Computation,Ncomponents, Type, Deltat, Stop_time , Stop_crit);
X_0 = 0;
Y_0 = 0;
d = 4;
V_0 = 2;

Example_Timepotential = @(t,x,y) (ramp_potential(t));
Physics2D = TimePotential_Var2d(Method, Physics2D, Example_Timepotential );
Save_Solution = 1;
Outputs = OutputsINI_Var2d(Method,Save_Solution);
Printing = 1;
Evo = 10;
Draw = 1;
Print = Print_Var2d(Printing,Evo,Draw);
[Phi,Outputs]= GPELab2d(Phi_0,Method,Geometry2D,Physics2D,Outputs,[],Print);