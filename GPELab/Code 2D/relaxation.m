Computation = 'Dynamic';
Ncomponents = 1;
Type = 'Relaxation';
Deltat = 1e-3;
Stop_time = 1.5;
Stop_crit = [];
Method = Method_Var2d(Computation,Ncomponents, Type, Deltat, Stop_time ,Stop_crit);
save = 1;
Evo_save = 100;
Outputs = OutputsINI_Var2d(Method,save,Evo_save);
Printing = 1;
Evo = 15;
Draw = 1;
Print = Print_Var2d(Printing,Evo,Draw);

[Phi,Outputs]= GPELab2d(Phi_3,Method,Geometry2D,Physics2D,Outputs,[],Print);