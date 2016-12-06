Computation = 'Dynamic';
Ncomponents = 1;
Type = 'Splitting';
Deltat = 1e-2;
t_hold = 4;
Stop_time = t_hold;
Method = Method_Var2d(Computation,Ncomponents, Type, Deltat, Stop_time);
Method.Splitting = 'Lie';
xmin = -2;
xmax = 2;
ymin = -2;
ymax = 2;
Nx = 2^10+1;
Ny = 2^10+1;
Geometry2D = Geometry2D_Var2d(xmin,xmax, ymin,ymax, Nx, Ny);
Delta = 0.5;
Beta = 3000;



v = 50;

k = 10;

Physics2D = TimePotential_Var2d(Method,Physics2D,@(t,x,y)(v*sin(k*x).^2));

Physics2D = Potential_Var2d(Method, Physics2D,@(x,y)(1/2)*(x.^2+y.^2));

Physics2D = Nonlinearity_Var2d(Method, Physics2D,@(phi,x,y) abs(phi).^2);

Outputs_iterations = 100;
Outputs_save = 1;
Outputs = OutputsINI_Var2d(Method,Outputs_iterations,Outputs_save);
Printing = 1;
Evo = 200;
Draw = 1;
Print = Print_Var2d(Printing,Evo,Draw);
[Phi_3,Outputs]= GPELab2d(Phi_2,Method,Geometry2D,Physics2D,Outputs,[],Print);