%% Creation of the figure variable
%% INPUTS(OPTIONAL):
%%          map: Choice of the mapcolor (char)
%% OUTPUT:
%%          Figure: Structure containing variables concerning the figures (structure)

function [Figure] = Figure_Var2d(varargin)

%% Figure's inputs
valid_Map = {'jet','hsv','hot','cool','spring','summer','autumn','winter','gray','bone','copper','pink','lines'}; % List of valid inputs for the map type
Analyse_Var = inputParser; % Creating the parser
Analyse_Var.addOptional('map','jet',@(x)ischar(validatestring(x,valid_Map))); % Optional input 'map' with default value 'jet'
Analyse_Var.addOptional('axis',[]); % Optional input 'caxis' with default value []
Analyse_Var.parse(varargin{:}); % Analysing the inputs

%% Contructing the figure structure
Figure.map = Analyse_Var.Results.map; % Storing the 'map' input, which is the colormap
Figure.axis = Analyse_Var.Results.axis; % Storing the 'caxis' input, which is the caxis
Figure.x = 'x';
Figure.y = 'y';
Figure.z = 'z';