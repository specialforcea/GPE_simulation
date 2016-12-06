%% Draw wave function's slice
%% INPUTS:
%%          phi: Wave function (matrix)
%%          Method: Structure containing variables concerning the method (structure)
%%          Geometry3D: Structure containing variables concerning the geometry of the problem in 3D (structure) (see Geometry3D_Var3d.m)
%%          Figure: Structure containing variables concerning the figures (structure)

function draw_slice_3d(phi,Geometry3D,Figure)
figure(Figure.label); % Number of the figure
clf(Figure.label); % Clear the figure
slice(Geometry3D.X,Geometry3D.Y,Geometry3D.Z,phi,Figure.Sy,Figure.Sx,Figure.Sz); % Computing the slice
shading interp; % Setting shading
colormap(Figure.map); % Setting the colormap
colorbar; % Adding the colorbar
xlabel(Figure.x); % Setting x-axis label
ylabel(Figure.y); % Setting y-axis label
title(Figure.title); % Setting title of the figure
drawnow;