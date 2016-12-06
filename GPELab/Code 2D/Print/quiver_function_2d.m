%% Draw 2D vector field
%% INPUTS:
%%          Density: density (matrix)
%%          Vecx, Vecy: real matrices (matrix)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          Geometry2D: Structure containing variables concerning the geometry of the problem in 2D (structure) (see Geometry2D_Var2d.m)
%%          Figure: Structure containing variables concerning the figures (structure) (see Figure_Var2d.m)

function quiver_function_2d(Density,Vecx,Vecy,Geometry2D,Figure)
figure(Figure.label); % Setting the number of the figure
clf(Figure.label); % Clear figure
hold on
quiver(Geometry2D.X,Geometry2D.Y,Vecx,Vecy); % Drawing function
contour(Geometry2D.X,Geometry2D.Y,Density); % Drawing contour density
axis equal
axis tight
shading interp; % Setting shading
colormap(Figure.map); % Setting colormap
colorbar; % Setting colorbar
xlabel(Figure.x); % Setting x-axis label
ylabel(Figure.y); % Setting y-axis label
title(Figure.title); % Setting title of the figure
if (isvector(Figure.axis)==1)
    caxis(Figure.axis)
end
drawnow; % Drawing
hold off