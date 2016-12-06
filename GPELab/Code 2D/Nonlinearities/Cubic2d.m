function [CubicNonlinearity] = Cubic2d(Method)

CubicNonlinearity = cell(Method.Ncomponents);
for n = 1:Method.Ncomponents
        CubicNonlinearity{n,n} =  @(Phi,X,Y) abs(Phi{n}).^2;
end