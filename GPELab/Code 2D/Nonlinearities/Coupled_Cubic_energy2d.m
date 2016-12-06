function [CoupledCubicEnergy] = Coupled_Cubic_energy2d(Beta)

CoupledCubicEnergy = cell(2);
CoupledCubicEnergy{1,1} = @(Phi,X,Y)(1/2)*Beta(1,1)*abs(Phi{1}).^2+(1/2)*Beta(1,2)*abs(Phi{2}).^2;
CoupledCubicEnergy{2,2} = @(Phi,X,Y)(1/2)*Beta(2,2)*abs(Phi{2}).^2+(1/2)*Beta(2,1)*abs(Phi{1}).^2;
CoupledCubicEnergy{1,2} = @(Phi,X,Y) 0;
CoupledCubicEnergy{2,1} = @(Phi,X,Y) 0;