%% Draw spin vector of the wave functions Phi
%% INPUTS:
%%          Phi: Wave functions (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          Geometry2D: Structure containing variables concerning the geometry of the problem in 2D (structure) (see Geometry2D_Var2d.m)
%%          Figure: Structure containing variables concerning the figures (structure) (see Figure_Var2d.m)
%% FUNCTIONS USED:
%%          draw_function_2d: To draw the wave function's square modulus and angle (line 19 and 23)

function Draw_Spin2d(Phi, Method, Geometry2D, Figure)

%% Computing the spin vector
Sx = real(conj(Phi{1}).*Phi{2})./(sqrt(abs(Phi{1}).^2+abs(Phi{2}).^2)); % spin vector x-component
Sy = imag(conj(Phi{1}).*Phi{2})./(sqrt(abs(Phi{1}).^2+abs(Phi{2}).^2)); % spin vector y-component
Density = sqrt(abs(Phi{1}).^2 + abs(Phi{2}).^2);

%% Printing the spin vector in the x-y plane
Figure.label = 1; % Number of the figure
Figure.title = 'Spin vector in the x-y plane (with density contour)'; % Storing title of the figure
quiver_function_2d(Density,Sx,Sy,Geometry2D,Figure); % Drawing the spin vector
