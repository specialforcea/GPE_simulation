    function InvA = Invc(A)
    %% Initializing zeros
    Zeros = zeros(size(A{1}));
        %% Computation of the inverse cell 
        if (length(A) ~= 1)
            InvA{1,1} = 1./A{1,1};
            for n = 1:length(A)-1
                %% Storing submatrices
                for k = 1:n
                    B{k} = A{k,n+1};
                    C{k} = A{n+1,k};
                end
                %% Computing intermediary matrices
                E = A{n+1,n+1};
                for k = 1:n
                    G{k} = Zeros;
                end
                for k = 1:n
                    for l = 1:n
                        G{k} = G{k} + C{l}.*InvA{l,k};
                    end
                end
                for k = 1:n
                    E = E - G{k}.*B{k};
                end
                E = 1./E;
                for k = 1:n
                    F{k} = Zeros;
                end
                for k = 1:n
                    for l = 1:n
                        F{k} = F{k} + InvA{k,l}.*B{l};
                    end
                end
                for k = 1:n
                    for l = 1:n
                        H{k,l} = E.*F{k}.*G{l};
                    end
                end
                %% Computing the intermediary inverse
                InvA{n+1,n+1} = E;
                for k = 1:n
                    InvA{k,n+1} = -E.*F{k};
                    InvA{n+1,k} = -E.*G{k};
                    for l = 1:n
                        InvA{k,l} = InvA{k,l}+H{k,l};
                    end
                end
                
            end
            %ELSEIF there is only one component
        elseif (length(A) == 1)
            InvA{1} = 1./A{1,1}; % Storing the full Thomas Fermi preconditioner
        end
    end
