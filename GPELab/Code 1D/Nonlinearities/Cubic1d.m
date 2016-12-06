function [CubicNonlinearity] = Cubic1d(Method)

CubicNonlinearity = cell(Method.Ncomponents);
for n = 1:Method.Ncomponents
        CubicNonlinearity{n,n} =  @(Phi,X) abs(Phi{n}).^2;
end