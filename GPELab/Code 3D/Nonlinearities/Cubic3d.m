function [CubicNonlinearity] = Cubic3d(Method)

CubicNonlinearity = cell(Method.Ncomponents);
for n = 1:Method.Ncomponents
        CubicNonlinearity{n,n} =  @(Phi,X,Y,Z) abs(Phi{n}).^2;
end