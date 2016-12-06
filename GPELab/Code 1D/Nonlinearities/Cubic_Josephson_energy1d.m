function [NLE] = Cubic_Josephson_energy1d(Beta11,Beta12,Beta22)

NLE = cell(2);
NLE{1,1} =  @(Phi,X) (1/2)*Beta11*abs(Phi{1}).^2+(1/2)*Beta12*abs(Phi{2}).^2;
NLE{2,2} =  @(Phi,X) (1/2)*Beta22*abs(Phi{2}).^2+(1/2)*Beta12*abs(Phi{1}).^2;
NLE{2,1} =  @(Phi,X) 0;
NLE{1,2} =  @(Phi,X) 0;