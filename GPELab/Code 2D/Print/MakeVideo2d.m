%% Creating a movie from the saved functions of a simulation
%% INPUTS:
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          Geometry2D: Structure containing variables concerning the geometry of the problem in 2D (structure) (see Geometry2D_Var2d.m)
%%          Outputs: Different outputs computed during the computation of the ground states (structure) (see OutputsINI_Var2D.m and FFTOutputs_Var2d.m)
%% INPUTS(OPTIONAL):
%%          Function: Function depending on (Phi,X,Y) which will be depicted in the video (cell array of functions or function)
%%          VideoName: Name of the video(s) (cell array of string or string)
%%          Figure: Structure containing variables concerning the figures (structure)

function MakeVideo2d(Method,Geometry2D,Outputs,varargin)
%% Setting default Function input
Default_Function = cell(1,Method.Ncomponents);
% FOR each components
for n = 1:Method.Ncomponents
    Default_Function{n} = @(phi,X,Y) abs(phi).^2; 
end
%% Analysis of inputs
Analyse_Var = inputParser; % Creating the parser
Analyse_Var.addOptional('Function',Default_Function); % Optional input 'Function'
Analyse_Var.addOptional('VideoName','MyVideo',@(x)(ischar(x)+iscell(x))); % Optional input 'Name'
Analyse_Var.addOptional('Figure', Figure_Var2d,@(x)isstruct(x)); % Optional input 'Figure' which must be a structure

%% Parsing inputs
% Parsing inputs
Analyse_Var.parse(varargin{:}); % Analysing the inputs
% Setting inputs
Function = Analyse_Var.Results.Function; %Storing the 'Function' input
Name = Analyse_Var.Results.VideoName; %Storing the 'Name' input
Figure = Analyse_Var.Results.Figure; %Storing the 'Name' input

% IF one has chosen to make a video for each component
if (iscell(Function) == 1)
    
for n = 1:Method.Ncomponents
% Setting figure and properties
figure(1);
clf(1);
pcolor(Geometry2D.X,Geometry2D.Y,Function{n}(Outputs.Solution{1}{n},Geometry2D.X,Geometry2D.Y))
axis equal
axis tight
shading interp; % Setting shading
colormap(Figure.map); % Setting colormap
colormap_tmp = caxis; % Store colormap setting
colorbar; % Setting colorbar
view(2); % Setting view
drawnow;
fprintf('Set the figure and then press any key \n')
pause
set(gcf, 'nextplot','replacechildren', 'Visible','off');
% Create a VideoObject
if (isempty(Name) == 0) && (iscell(Name) == 0)
    vidObj = VideoWriter(strcat(Name,num2str(n),'.avi'));
elseif (isempty(Name) == 0) && (iscell(Name) == 1)
    vidObj = VideoWriter(strcat(Name{n},'.avi'));
elseif (isempty(Name) == 1)
    vidObj = VideoWriter(strcat('MyVideo',num2str(n),'.avi'));
end
vidObj.Quality = 100;
vidObj.FrameRate = 35;
open(vidObj);

for m = 1:Outputs.Iterations
   pcolor(Geometry2D.X,Geometry2D.Y,Function{n}(Outputs.Solution{m}{n},Geometry2D.X,Geometry2D.Y))
   axis equal
   axis tight
   caxis(colormap_tmp)
   shading interp;
   view(2)
   writeVideo(vidObj, getframe(gcf));
end
close(vidObj);
end

% ELSEIF one has chosen to make a single video
elseif (iscell(Function) == 0)
% Setting figure and properties
figure(1);
clf(1);
surf(Geometry2D.X,Geometry2D.Y,Function(Outputs.Solution{1}{1},Geometry2D.X,Geometry2D.Y),'EdgeColor','none')
shading interp;
view(2)
drawnow;
fprintf('Set the figure and then press any key \n')
pause
set(gcf, 'nextplot','replacechildren', 'Visible','off');
% Create a VideoObject
if (iscell(Name) == 0)
    vidObj = VideoWriter(strcat(Name,'.avi'));
elseif (iscell(Name) == 1)
    vidObj = VideoWriter('MyVideo.avi');
end
vidObj.Quality = 100;
vidObj.FrameRate = 12;
open(vidObj);

for m = 1:Outputs.Iterations
   surf(Geometry2D.X,Geometry2D.Y,Function(Outputs.Solution{m},Geometry2D.X,Geometry2D.Y),'EdgeColor','none')
   shading interp;
   view(2)
   writeVideo(vidObj, getframe(gcf));
end
close(vidObj);
end