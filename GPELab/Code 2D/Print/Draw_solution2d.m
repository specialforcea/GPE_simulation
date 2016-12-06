%% Draw wave functions' square of the modulus and angle
%% INPUTS:
%%          Phi: Wave functions (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          Geometry2D: Structure containing variables concerning the geometry of the problem in 2D (structure) (see Geometry2D_Var2d.m)
%%          Figure: Structure containing variables concerning the figures (structure) (see Figure_Var2d.m)
%% FUNCTIONS USED:
%%          draw_function_2d: To draw the wave function's square modulus and angle (line 19 and 23)

function Draw_solution2d(Phi, Method, Geometry2D, Figure, varargin)
%% Setting default Function input
Default_Function = @(phi,X,Y) abs(phi).^2;
%% Analysis of inputs
Analyse_Var = inputParser; % Creating the parser
Analyse_Var.addOptional('Function',Default_Function); % Optional input 'Function'
Analyse_Var.addOptional('Function_Name','|\phi(x,y)|^2',@(x) (ischar(x)) || (iscell(x))); % Optional input 'Function_Name'

%% Parsing inputs
% Parsing inputs
Analyse_Var.parse(varargin{:}); % Analysing the inputs
% Setting inputs
Function = Analyse_Var.Results.Function; %Storing the 'Function' input
Function_Name = Analyse_Var.Results.Function_Name; %Storing the 'Function_Name' input

if (iscell(Function) == 0)
% FOR each component
for n = 1:Method.Ncomponents
    %% Printing the figure of the square of the modulus of wave function
    Figure.label = n; % Number of the figure
    Figure.title = strcat(Function_Name,' of component ', 32, num2str(n)); % Storing title of the figure
    draw_function_2d(Function(Phi{n},Geometry2D.X,Geometry2D.Y),Geometry2D,Figure); % Drawing the square of the modulus of the wave function
    %% Printing the figure of the angle of wave function
%     Figure.label = n+Method.Ncomponents; % Number of the figure
%     Figure.title = strcat('angle(\phi',32,'(x,y)) of component ', 32 ,num2str(n)); % Storing title of the figure
%     draw_function_2d(angle(Phi{n}),Geometry2D,Figure); % Drawing the angle of the wave function
end
elseif (iscell(Function) == 1)
% FOR each component
for n = 1:length(Function)
    %% Printing the figure of the square of the modulus of wave function
    Figure.label = n; % Number of the figure
    Figure.title = strcat(Function_Name{n}); % Storing title of the figure
    draw_function_2d(Function{n}(Phi,Geometry2D.X,Geometry2D.Y),Geometry2D,Figure); % Drawing the square of the modulus of the wave function
end
end