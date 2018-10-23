function M = interpolation_matrix(Nx,Ns)
factor = int32(Nx/Ns)+1;
M = zeros(Nx,Ns);
j = 1;
i = int32(1);
while i+factor < Nx
    M(i:i+factor-1, j)= ones(factor,1);
    i = i+factor;
    j = j + 1;
end

M(i:Nx,j) = ones(Nx-i+1,1);

end
    

