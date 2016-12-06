function [NL] = Cubic_Josephson1d(Beta11,Beta12,Beta22)

NL = cell(2);
NL{1,1} =  @(Phi,X) Beta11*abs(Phi{1}).^2+Beta12*abs(Phi{2}).^2;
NL{2,2} =  @(Phi,X) Beta22*abs(Phi{2}).^2+Beta12*abs(Phi{1}).^2;
NL{2,1} =  @(Phi,X) 0;
NL{1,2} =  @(Phi,X) 0;