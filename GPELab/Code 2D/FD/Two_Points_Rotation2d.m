%% Computation of the rotational operator with a two points finite difference scheme and Dirichlet conditions
%% INPUT:
%%          FDGeometry2D: Structure containing variables concerning the geometry of the problem in 2D in the FD context (structure) (see FDGeometry2D_Var2d.m)
%% OUTPUT:
%%          L_z: Discretized rotational with the FD scheme (matrix)

function [ L_z ] = Two_Points_Rotation2d(FDGeometry2D)
N2 = FDGeometry2D.Nx*FDGeometry2D.Ny; %Total number of internal grid points
%% Computation of the operator y\grad_x
yDiag = repmat(FDGeometry2D.Y(:,1),FDGeometry2D.Nx,1); %Computing the vector corresponding to the diagonal of the y\grad_x operator
yDiagx = spdiags([yDiag],FDGeometry2D.Ny,N2,N2); %Computing the diagonal of the y\grad_x operator
yGx = (1/(2*FDGeometry2D.dx))*(yDiagx-yDiagx'); %Assembling the y\grad_x operator

%% Computation of the operator x\grad_y
xDiagy = [];
for i=1:FDGeometry2D.Nx
xDiagy = [xDiagy;[0;ones(FDGeometry2D.Ny-1,1)]*FDGeometry2D.X(1,i)]; %Computing the vector corresponding to the diagonal of the x\grad_y operator
end
xDiagy = spdiags(xDiagy,1,N2,N2); %Computing the diagonal of the x\grad_y operator
xGy = (1/(2*FDGeometry2D.dy))*(xDiagy-xDiagy'); %Assembling the x\grad_y operator

%% Computation of the rotational operator
L_z = -1i*(xGy-yGx); %Assembling the rotational operator