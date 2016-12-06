function AB = Mulc(A,B)
DimA = size(A);
DimB = size(B);
dimA = size(A{1});

if iscell(B)

for n = 1:DimA(1)
    for m = 1:DimB(2)
        AB{n,m} = zeros(dimA);
    end
end

for n = 1:DimA(1)
    for m = 1:DimB(2)
        for l = 1:DimA(2)
            AB{n,m} = AB{n,m} + A{n,l}.*B{l,m};
        end
    end
end

elseif isscalar(B)
    
for n = 1:DimA(1)
    for m = 1:DimA(2)
        AB{n,m} = B*A{n,m};
    end
end

end