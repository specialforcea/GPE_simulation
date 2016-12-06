function [CubicNonlinearity] = Cubic_energy3d(Method)

CubicNonlinearity = cell(Method.Ncomponents);
for n = 1:Method.Ncomponents
        CubicNonlinearity{n,n} =  @(Phi,X,Y,Z) (1/2)*abs(Phi{n}).^2;
end