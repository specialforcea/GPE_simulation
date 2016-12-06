%% Draw wave functions' square of the modulus and angle in time
%% INPUTS:
%%          Outputs: Different outputs computed during the computation of the ground states (structure) (see OutputsINI_Var1D.m and FFTOutputs_Var1d.m)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var1d.m)
%%          Geometry1D: Structure containing variables concerning the geometry of the problem in 1D (structure) (see Geometry1D_Var1d.m)
%%          Figure: Structure containing variables concerning the figures (structure) (see Figure_Var1d.m)

function Draw_Timecontour1d(Outputs, Method, Geometry1D, Figure)
% FOR each component
for n = 1:Method.Ncomponents
    %% Building matrix associated with the evolution of the wave function
    PHI = zeros(Geometry1D.Nx,Outputs.Iterations);
    for m = 1:Outputs.Iterations
        PHI(:,m) = Outputs.Solution{m}{n};
    end
    Time = 0 : Method.Deltat*Outputs.Evo_outputs : Method.Stop_time;
    %% Printing the figure of the square of the modulus of wave function
    Figure.label = n; % Number of the figure
    Figure.title = strcat('|\phi(x)',32,'|^2 of component ', 32, num2str(n)); % Storing title of the figure
    figure(Figure.label); % Setting the number of the figure
    clf(Figure.label); % Clear figure
    V = [0.01,0.1,0.5,1];
    contour(Time,Geometry1D.X,abs(PHI).^2,V); % Drawing function
    colormap('jet'); % Setting colormap
    colorbar; % Setting colorbar
    ylabel(Figure.x); % Setting x-axis label
    xlabel('Time')
    title(Figure.title); % Setting title of the figure
    view(2)
    drawnow; % Drawing
end