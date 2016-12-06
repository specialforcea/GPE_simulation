%% Creation of the 1D geometry structure
%% INPUT:
%%          xmin,xmax: Lengths from the origin of the square domain in the x direction (double)
%%          Nx: Number of points for the discretization of the domain in the x direction (double)
%% OUTPUT:
%%          Geometry1D: Structure containing variables concerning the 1D geometry (structure)

function [Geometry1D] = Geometry1D_Var1d(varargin)
%% Analysis of the inputs
Analyse_Var = inputParser; % Creating the parser
Analyse_Var.addOptional('xmin',-10,@(x)isscalar(x)); % Optional input 'xmin' with default value -10
Analyse_Var.addOptional('xmax',10,@(x)isscalar(x)); % Optional input 'xmax' with default value 10
Analyse_Var.addOptional('Nx',2^7+1,@(x)isposintscalar(x)); % Optional input 'Nx' with default value 2e7

%% Parsing inputs and creating the Method structure
% Parsing inputs
Analyse_Var.parse(varargin{:}); % Analysing the inputs
% Contructing the Method structure
xmin = Analyse_Var.Results.xmin; %Storing the 'xmin' input
xmax = Analyse_Var.Results.xmax; %Storing the 'xmax' input
Nx = ceil(Analyse_Var.Results.Nx); %Storing the 'Nx' input

%% Computation of the geometric variables
Geometry1D.Lx = xmax-xmin; % Size of the domain in the x direction
Geometry1D.Nx = Nx; % Number of points in the x direction
Geometry1D.dx = (xmax-xmin)/(Nx-1); % Spatial step along x
Geometry1D.X = (xmin:Geometry1D.dx:xmax)'; % Coordinates of the points in the x,y,z direction