function y = JSD(p,q,Nx,deltax)
M = (p +q)./2;
y = 1/2*(integr(p.*log2(p./M),Nx,deltax) + integr(q.*log2(q./M),Nx,deltax));
end