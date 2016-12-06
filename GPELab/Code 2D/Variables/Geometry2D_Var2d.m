%% Creation of the 2D geometry structure
%% INPUT:
%%          xmin,xmax: Lengths from the origin of the square domain in the x direction (double)
%%          ymin,ymax: Lengths from the origin of the square domain in the y direction (double)
%%          Nx,Ny: Number of points for the discretization of the domain in the x,y direction (double)
%% OUTPUT:
%%          Geometry2D: Structure containing variables concerning the 2D geometry (structure)

function [Geometry2D] = Geometry2D_Var2d(varargin)
%% Analysis of the inputs
Analyse_Var = inputParser; % Creating the parser
Analyse_Var.addOptional('xmin',-10,@(x)isscalar(x)); % Optional input 'xmin' with default value -10
Analyse_Var.addOptional('xmax',10,@(x)isscalar(x)); % Optional input 'xmax' with default value 10
Analyse_Var.addOptional('ymin',-10,@(x)isscalar(x)); % Optional input 'ymin' with default value -10
Analyse_Var.addOptional('ymax',10,@(x)isscalar(x)); % Optional input 'ymax' with default value 10
Analyse_Var.addOptional('Nx',2^7+1,@(x)isposintscalar(x)); % Optional input 'Nx' with default value 2^7+1
Analyse_Var.addOptional('Ny',2^7+1,@(x)isposintscalar(x)); % Optional input 'Ny' with default value 2^7+1

%% Parsing inputs and creating the Method structure
% Parsing inputs
Analyse_Var.parse(varargin{:}); % Analysing the inputs
% Contructing the Method structure
xmin = Analyse_Var.Results.xmin; %Storing the 'xmin' input
xmax = Analyse_Var.Results.xmax; %Storing the 'xmax' input
ymin = Analyse_Var.Results.ymin; %Storing the 'ymin' input
ymax = Analyse_Var.Results.ymax; %Storing the 'ymax' input
Nx = ceil(Analyse_Var.Results.Nx); %Storing the 'Nx' input
Ny = ceil(Analyse_Var.Results.Ny); %Storing the 'Ny' input

%% Computation of the geometric variables
Geometry2D.Lx = xmax-xmin; % Size of the domain in the x direction
Geometry2D.Ly = ymax-ymin; % Size of the domain in the y direction
Geometry2D.Nx = Nx; % Number of points in the x direction
Geometry2D.Ny = Ny; % Number of points in the y direction
Geometry2D.N2 = Nx*Ny; % Total number of points in the domain
Geometry2D.dx = (xmax-xmin)/(Nx-1); % Spatial step along x
Geometry2D.dy = (ymax-ymin)/(Ny-1); % Spatial step along y
[Geometry2D.X,Geometry2D.Y] = meshgrid(xmin:Geometry2D.dx:xmax,ymin:Geometry2D.dy:ymax); % Coordinates of the points in the x,y,z direction