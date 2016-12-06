%% Draw wave function
%% INPUTS:
%%          phi: Wave function (matrix)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var1d.m)
%%          Geometry1D: Structure containing variables concerning the geometry of the problem in 1D (structure) (see Geometry1D_Var1d.m)
%%          Figure: Structure containing variables concerning the figures (structure) (see Figure_Var1d.m)

function draw_function_1d(phi,Geometry1D,Figure)
figure(Figure.label); % Setting the number of the figure
clf(Figure.label); % Clear figure
plot(Geometry1D.X,phi); % Drawing function
xlabel(Figure.x); % Setting x-axis label
title(Figure.title); % Setting title of the figure
axis([-15,15,0,1e-1])
drawnow; % Drawing