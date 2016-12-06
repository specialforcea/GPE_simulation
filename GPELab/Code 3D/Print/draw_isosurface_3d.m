%% Draw wave function's isosurface
%% INPUTS:
%%          phi: Wave function (matrix)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var3d.m)
%%          Geometry3D: Structure containing variables concerning the geometry of the problem in 3D (structure) (see Geometry3D_Var3d.m)
%%          Figure: Structure containing variables concerning the figures (structure) (see Figure_Var3d.m)

function draw_isosurface_3d(phi,Geometry3D,Figure)
figure(Figure.label); % Number of the figure
clf(Figure.label); % Clear the figure
hpatch = patch(isosurface(Geometry3D.X,Geometry3D.Y,Geometry3D.Z,phi,Figure.iso)); % Computing the patch
isonormals(Geometry3D.X,Geometry3D.Y,Geometry3D.Z,phi,hpatch); % Computing the isosnormal
set(hpatch,'facecolor',Figure.color,'edgecolor','none'); % Setting the color of faces and edges
daspect(Figure.aspect); % Setting the data aspect ratio
view(Figure.view); % Setting the angle of view
axis tight; % Setting the axis limits to the range of the data
grid on; % Adding major grid lines
camlight; % Setting light
axis  equal; % Setting the aspect ratio so that increments on the x,y,z axis are of equal size
title(Figure.title) % Adding title
xlabel(Figure.x); % Adding text beside the x-axis
xlabel(Figure.y); % Adding text beside the y-axis
xlabel(Figure.z); % Adding text beside the z-axis
lighting gouraud; % Setting the lighting to gouraud
alpha(Figure.alpha); % Setting transparency
drawnow;