function [CoupledCubicNonlinearity] = Coupled_Cubic3d(Beta)

CoupledCubicNonlinearity = cell(2);
CoupledCubicNonlinearity{1,1} = @(Phi,X,Y,Z)Beta(1,1)*abs(Phi{1}).^2 + Beta(1,2)*abs(Phi{2}).^2;
CoupledCubicNonlinearity{2,2} = @(Phi,X,Y,Z)Beta(2,2)*abs(Phi{2}).^2 + Beta(2,1)*abs(Phi{1}).^2;
CoupledCubicNonlinearity{1,2} = @(Phi,X,Y,Z) 0;
CoupledCubicNonlinearity{2,1} = @(Phi,X,Y,Z) 0;