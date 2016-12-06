function AB = Addc(A,B)

DimA = size(A);

if iscell(B)
    
for n = 1:DimA(1)
    for m = 1:DimA(2)
        AB{n,m} = A{n,m} + B{n,m};
    end
end

elseif isscalar(B)
    
for n = 1:DimA(1)
    for m = 1:DimA(2)
        AB{n,m} = A{n,m} + B;
    end
end

end