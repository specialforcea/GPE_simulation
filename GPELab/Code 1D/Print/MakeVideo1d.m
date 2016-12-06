%% Creating a movie from the saved functions of a simulation
%% INPUTS:
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var1d.m)
%%          Geometry1D: Structure containing variables concerning the geometry of the problem in 2D (structure) (see Geometry1D_Var1d.m)
%%          Outputs: Different outputs computed during the computation of the ground states (structure) (see OutputsINI_Var1D.m and FFTOutputs_Var1d.m)
%% INPUTS(OPTIONAL):
%%          MultiVideo: Must either be 1 if one wants to print each component on a single frame or 0 if not (integer)
%%          Function: Function depending on (Phi,X) which will be depicted in the video (cell array of functions or function)
%%          VideoName: Name of the video(s) (cell array of string or string)

function MakeVideo1d(Method,Geometry1D,Outputs,varargin)
%% Setting default Function input
Default_Function = cell(1,Method.Ncomponents);
% FOR each components
for n = 1:Method.Ncomponents
    Default_Function{n} = @(phi,X) abs(phi).^2; 
end

%% Setting colors
Colors = {'b','r','g','c','m','y','k'};

%% Analysis of inputs
Analyse_Var = inputParser; % Creating the parser
Analyse_Var.addOptional('MultiVideo', 0); % Optional input 'MultiVideo'
Analyse_Var.addOptional('Function',Default_Function); % Optional input 'Function'
Analyse_Var.addOptional('VideoName','MyVideo',@(x)(ischar(x)+iscell(x))); % Optional input 'Name'

%% Parsing inputs
% Parsing inputs
Analyse_Var.parse(varargin{:}); % Analysing the inputs
% Setting inputs
MultiVideo = Analyse_Var.Results.MultiVideo; %Storing the 'MultiVideo' input
Function = Analyse_Var.Results.Function; %Storing the 'Computation' input
Name = Analyse_Var.Results.VideoName; %Storing the 'Ncomponents' input

% IF one has chosen to print each component on a different frame
if (MultiVideo == 0)

% IF one has chosen to make a video for each component
if (iscell(Function) == 1)

for n = 1:Method.Ncomponents
% Setting figure and properties
figure(2);
clf(2);
plot(Geometry1D.X,Function{n}(Outputs.Solution{1}{n},Geometry1D.X))
drawnow;
fprintf('Set the figure and then press any key \n')
pause
% Create a VideoObject
if (isempty(Name) == 0) && (iscell(Name) == 0)
    vidObj = VideoWriter(strcat(Name,num2str(n),'.avi'));
elseif (isempty(Name) == 0) && (iscell(Name) == 1)
    vidObj = VideoWriter(strcat(Name{n},'.avi'));
elseif (isempty(Name) == 1)
    vidObj = VideoWriter(strcat('MyVideo',num2str(n),'.avi'));
end
vidObj.Quality = 100;
vidObj.FrameRate = 12;
open(vidObj);

for m = 1:Outputs.Iterations
   plot(Geometry1D.X,Function{n}(Outputs.Solution{m}{n},Geometry1D.X))
   writeVideo(vidObj, getframe);
end
close(vidObj);
end

% ELSEIF one has chosen to make a single video
elseif (iscell(Function) == 0)
% Setting figure and properties
figure(2);
clf(2);
plot(Geometry1D.X,Function(Outputs.Solution{1}{1},Geometry1D.X))
drawnow;
fprintf('Set the figure and then press any key \n')
pause
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
   plot(Geometry1D.X,Function(Outputs.Solution{m},Geometry1D.X))
   writeVideo(vidObj, getframe);
end
close(vidObj);
end

% ELSE IF one has chosen to print each component on a single frame
elseif (MultiVideo == 1)
    
% IF one has chosen to make a video for each component
if (iscell(Function) == 1)


% Setting figure and properties
figure(2);
clf(2);
hold on
for n = 1:Method.Ncomponents
plot(Geometry1D.X,Function{n}(Outputs.Solution{1}{n},Geometry1D.X),Colors{mod(n,7)})
drawnow;
end
fprintf('Set the figure and then press any key \n')
pause
% Create a VideoObject
if (isempty(Name) == 0) && (iscell(Name) == 0)
    vidObj = VideoWriter(strcat(Name,'.avi'));
elseif (isempty(Name) == 0) && (iscell(Name) == 1)
    vidObj = VideoWriter(strcat(Name{1},'.avi'));
elseif (isempty(Name) == 1)
    vidObj = VideoWriter(strcat('MyVideo','.avi'));
end
vidObj.Quality = 100;
vidObj.FrameRate = 12;
open(vidObj);

for m = 1:Outputs.Iterations
   for n = 1:Method.Ncomponents
   plot(Geometry1D.X,Function{n}(Outputs.Solution{m}{n},Geometry1D.X),Colors{mod(n,7)})
   hold on
   end
   hold off
   writeVideo(vidObj, getframe);
end
close(vidObj);
end

% ELSEIF one has chosen to make a single video
elseif (iscell(Function) == 0)
% Setting figure and properties
figure(2);
clf(2);
plot(Geometry1D.X,Function(Outputs.Solution{1}{1},Geometry1D.X))
drawnow;
fprintf('Set the figure and then press any key \n')
pause
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
   plot(Geometry1D.X,Function(Outputs.Solution{m},Geometry1D.X))
   writeVideo(vidObj, getframe);
end
close(vidObj);
end
    
end