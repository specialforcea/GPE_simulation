%% Creation of the figure variable
%% INPUTS:
%%          map: Choice of the mapcolor (char)
%%          alpha: Transparency of the isosurface (double)
%%          view: Set the angle of view in the 3D plot (double) (see view.m)
%%          iso: Set the isovalue for the isosurface (double)
%% OUTPUT:
%%          Figure: Structure containing variables concerning the figures (structure)

function [Figure] = Figure_Var3d(varargin)

%% Figure's inputs
valid_Map = {'jet','hsv','hot','cool','spring','summer','autumn','winter','gray','bone','copper','pink','lines'}; % List of valid inputs for the map type
Analyse_Var = inputParser; % Creating the parser
Analyse_Var.addOptional('view',3,@(x)isnumeric(x)); % Optional input 'view' with default value [85,10]
Analyse_Var.addOptional('iso',0.001,@(x)isposrealscalar(x)); % Optional input 'iso' with default value 0.001
Analyse_Var.addOptional('alpha',0.6,@(x)isposrealscalar(x)); % Optional input 'alpha' with default value 1/2
Analyse_Var.addOptional('aspect',[1,1,1],@(x)isnumeric(x)); % Optional input 'aspect' with default value [1,1,1]
Analyse_Var.addOptional('Sx', [0],@(x)isnumeric(x)); % Optional input 'Sx' with default value [0,0,0]
Analyse_Var.addOptional('Sy', [0],@(x)isnumeric(x)); % Optional input 'Sy' with default value [0,0,0]
Analyse_Var.addOptional('Sz', [0],@(x)isnumeric(x)); % Optional input 'Sz' with default value [0,0,0]
Analyse_Var.addOptional('map','jet',@(x)ischar(validatestring(x,valid_Map))); % Optional input 'map' with default value 'jet'
Analyse_Var.parse(varargin{:}); % Analysing the inputs

%% Contructing the figure structure
Figure.map = Analyse_Var.Results.map; % Storing the 'map' input, which is the colormap
Figure.alpha = Analyse_Var.Results.alpha; % Storing the 'alpha' input, which is the transparency
Figure.aspect = Analyse_Var.Results.aspect; % Storing the 'aspect' input, which is the data ratio aspect
Figure.view = Analyse_Var.Results.view; % Storing the 'view' input, which is the angle of view
Figure.iso = Analyse_Var.Results.iso; % Storing the 'iso' input, which is the isovalue for the isosurface
Figure.Sx = Analyse_Var.Results.Sx; % Storing the 'Sx' input, which is the coordinates of the points along the x axis where to draw the slices
Figure.Sy = Analyse_Var.Results.Sy; % Storing the 'Sy' input, which is the coordinates of the points along the y axis where to draw the slices
Figure.Sz = Analyse_Var.Results.Sz; % Storing the 'Sz' input, which is the coordinates of the points along the z axis where to draw the slices

Figure.x = 'x';
Figure.y = 'y';
Figure.z = 'z';

Figure.color = 'blue';