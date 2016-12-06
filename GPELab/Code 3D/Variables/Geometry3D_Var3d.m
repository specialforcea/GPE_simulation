%% Creation of the 3D geometry structure
%% INPUT:
%%          xmin,xmax: Lengths from the origin of the square domain in the x direction (double)
%%          ymin,ymax: Lengths from the origin of the square domain in the y direction (double)
%%          zmin,zmax: Lengths from the origin of the square domain in the z direction (double)
%%          Nx,Ny;Nz: Number of points for the discretization of the domain in the x,y,z direction (double)
%% OUTPUT:
%%          Geometry3D: Structure containing variables concerning the 3D geometry (structure)

function [Geometry3D] = Geometry3D_Var3d(varargin)
%% Analysis of the inputs
Analyse_Var = inputParser; % Creating the parser
Analyse_Var.addOptional('xmin',-10,@(x)isscalar(x)); % Optional input 'xmin' with default value -10
Analyse_Var.addOptional('xmax',10,@(x)isscalar(x)); % Optional input 'xmax' with default value 10
Analyse_Var.addOptional('ymin',-10,@(x)isscalar(x)); % Optional input 'ymin' with default value -10
Analyse_Var.addOptional('ymax',10,@(x)isscalar(x)); % Optional input 'ymax' with default value 10
Analyse_Var.addOptional('zmin',-10,@(x)isscalar(x)); % Optional input 'zmin' with default value -10
Analyse_Var.addOptional('zmax',10,@(x)isscalar(x)); % Optional input 'zmax' with default value 10
Analyse_Var.addOptional('Nx',2e7+1,@(x)isposintscalar(x)); % Optional input 'Nx' with default value 2e7
Analyse_Var.addOptional('Ny',2e7+1,@(x)isposintscalar(x)); % Optional input 'Ny' with default value 2e7
Analyse_Var.addOptional('Nz',2e7+1,@(x)isposintscalar(x)); % Optional input 'Nz' with default value 2e7

%% Parsing inputs and creating the Method structure
% Parsing inputs
Analyse_Var.parse(varargin{:}); % Analysing the inputs
% Contructing the Method structure
xmin = Analyse_Var.Results.xmin; %Storing the 'xmin' input
xmax = Analyse_Var.Results.xmax; %Storing the 'xmax' input
ymin = Analyse_Var.Results.ymin; %Storing the 'ymin' input
ymax = Analyse_Var.Results.ymax; %Storing the 'ymax' input
zmin = Analyse_Var.Results.zmin; %Storing the 'zmin' input
zmax = Analyse_Var.Results.zmax; %Storing the 'zmax' input
Nx = ceil(Analyse_Var.Results.Nx); %Storing the 'Nx' input
Ny = ceil(Analyse_Var.Results.Ny); %Storing the 'Ny' input
Nz = ceil(Analyse_Var.Results.Nz); %Storing the 'Ny' input

%% Computation of the geometric variables
Geometry3D.Lx = xmax-xmin; % Size of the domain in the x direction
Geometry3D.Ly = ymax-ymin; % Size of the domain in the y direction
Geometry3D.Lz = zmax-zmin; % Size of the domain in the z direction
Geometry3D.Nx = Nx; % Number of points in the x direction
Geometry3D.Ny = Ny; % Number of points in the y direction
Geometry3D.Nz = Nz; % Number of points in the z direction
Geometry3D.N3 = Nx*Ny*Nz; % Total number of points in the domain
Geometry3D.dx = (xmax-xmin)/(Nx-1); % Spatial step along x
Geometry3D.dy = (ymax-ymin)/(Ny-1); % Spatial step along y
Geometry3D.dz = (zmax-zmin)/(Nz-1); % Spatial step along z
[Geometry3D.X,Geometry3D.Y,Geometry3D.Z] = meshgrid(xmin:Geometry3D.dx:xmax,ymin:Geometry3D.dy:ymax,zmin:Geometry3D.dz:zmax); % Coordinates of the points in the x,y,z direction