%% Draw wave function's projection on the x,y plane
%% INPUTS:
%%          phi: Wave function (matrix)
%%          Method: Structure containing variables concerning the method (structure)
%%          Geometry3D: Structure containing variables concerning the geometry of the problem in 3D (structure) (see Geometry3D_Var3d.m)
%%          Figure: Structure containing variables concerning the figures (structure)

function Draw_projection_3d(Phi,Method,Geometry3D,Figure)

for n = 1:Method.Ncomponents
    %% Printing the figure of the square of the modulus of wave function
    figure(n); % Number of the figure
    clf(n); % Clear the figure
    Phi_projz{n} = sqrt(Geometry3D.dz)*sqrt(sum(abs(Phi{n}).^2,3));
    surf(Geometry3D.X(:,:,1),Geometry3D.Y(:,:,1),Phi_projz{n},'EdgeColor','none'); % Drawing function
    shading interp; % Setting shading
    colormap(Figure.map); % Setting colormap
    colorbar; % Setting colorbar
    view(2); % Setting view
    xlabel(Figure.x); % Setting x-axis label
    ylabel(Figure.y); % Setting y-axis label
    title(strcat('|\phi(x,y,.)|^2 of component ', 32, num2str(n))); % Setting title of the figure
    drawnow; % Drawing
end


