function F = expc(A)
%EXPC   Exponential cell.
%   expc(X) is the cell exponential of X.  expc is computed using
%   a scaling and squaring algorithm with a Pade approximation.

%   Reference:
%   N. J. Higham, The scaling and squaring method for the matrix
%   exponential revisited. SIAM J. Matrix Anal. Appl.,
%   26(4) (2005), pp. 1179-1193.

[m_vals, theta, classA] = expmchk; % Initialization
DimA = size(A);
for k = 1:DimA(1)
    for l = 1:DimA(2)
        MaxA(k,l) = max(max(max(A{k,l})));
    end
end

normA = norm(MaxA,1);

if normA <= theta(end)
    % no scaling and squaring is required.
    for i = 1:length(m_vals)
        if normA <= theta(i)
            F = PadeApproximantOfDegree(m_vals(i));
            break;
        end
    end
else

    [t, s] = log2(normA/theta(end));
    s = s - (t == 0.5); % adjust s if normA/theta(end) is a power of 2.
    A = Mulc(A,1/2^s);    % Scaling
    F = PadeApproximantOfDegree(m_vals(end));
    for i = 1:s
        F = Mulc(F,F);  % Squaring
    end
end
% End of expc

%%%%Nested Functions%%%%
    function [m_vals, theta, classA] = expmchk
        %EXPMCHK Check the class of input A and
        %    initialize M_VALS and THETA accordingly.
        classA = class(A{1});
        switch classA
            case 'double'
                m_vals = [3 5 7 9 13];
                % theta_m for m=1:13.
                theta = [%3.650024139523051e-008
                         %5.317232856892575e-004
                          1.495585217958292e-002  % m_vals = 3
                         %8.536352760102745e-002
                          2.539398330063230e-001  % m_vals = 5
                         %5.414660951208968e-001
                          9.504178996162932e-001  % m_vals = 7
                         %1.473163964234804e+000
                          2.097847961257068e+000  % m_vals = 9
                         %2.811644121620263e+000
                         %3.602330066265032e+000
                         %4.458935413036850e+000
                          5.371920351148152e+000];% m_vals = 13
            case 'single'
                m_vals = [3 5 7];
                % theta_m for m=1:7.
                theta = [%8.457278879935396e-004
                         %8.093024012430565e-002
                          4.258730016922831e-001  % m_vals = 3
                         %1.049003250386875e+000
                          1.880152677804762e+000  % m_vals = 5
                         %2.854332750593825e+000
                          3.925724783138660e+000];% m_vals = 7
            otherwise
                error(message('MATLAB:expm:inputType'))
        end
    end

    function F = PadeApproximantOfDegree(m)
        %PADEAPPROXIMANTOFDEGREE  Pade approximant to exponential.
        %   F = PADEAPPROXIMANTOFDEGREE(M) is the degree M diagonal
        %   Pade approximant to EXP(A), where M = 3, 5, 7, 9 or 13.
        %   Series are evaluated in decreasing order of powers, which is
        %   in approx. increasing order of maximum norms of the terms.
        n = length(A);
        c = getPadeCoefficients;
        Eyec = cell(DimA);
        Zeroc = cell(DimA);
        for k = 1:DimA(1)
            for l = 1:DimA(2)
                if (k == l)
                    Eyec{k,l} = 1;
                else
                    Eyec{k,l} = 0;
                end
                Zeroc{k,l} = 0;
            end
        end
        % Evaluate Pade approximant.
        switch m

            case {3, 5, 7, 9}

                Apowers = cell(ceil((m+1)/2),1);
                Apowers{1} = Eyec;
                Apowers{2} = Mulc(A,A);
                for j = 3:ceil((m+1)/2)
                    Apowers{j} = Mulc(Apowers{j-1},Apowers{2});
                end
                U = Zeroc; V = Zeroc;

                for j = m+1:-2:2
                    U = Addc(U,Mulc(Apowers{j/2},c(j)));
                end
                U = Mulc(A,U);
                for j = m:-2:1
                    V = Addc(V,Mulc(Apowers{(j+1)/2},c(j)));
                end
                F = Mulc(Invc(Addc(Mulc(U,-1),V)),Addc(U,V));

            case 13

                % For optimal evaluation need different formula for m >= 12.
                A2 = Mulc(A,A); A4 = Mulc(A2,A2); A6 = Mulc(A2,A4);
                U1 = Addc(Addc(Mulc(A6,c(14)), Mulc(A4,c(12))), Mulc(A2,c(10)));
                U2 = Addc(Addc(Addc(Mulc(A6,c(8)), Mulc(A4,c(6))), Mulc(A2,c(4))), Mulc(Eyec,c(2)));
                U = Mulc(A, (Addc(Mulc(A6,U1),U2)));
                V1 = Addc(Addc(Mulc(A6,c(13)), Mulc(A4,c(11))), Mulc(A2,c(9)));
                V2 = Addc(Addc(Addc(Mulc(A6,c(7)), Mulc(A4,c(5))), Mulc(A2,c(3))), Mulc(Eyec,c(1)));
                V = Addc(Mulc(A6,V1), V2);
                F = Mulc(Invc(Addc(Mulc(U,-1),V)),Addc(U,V));

        end

        function c = getPadeCoefficients
            % GETPADECOEFFICIENTS Coefficients of numerator P of Pade approximant
            %    C = GETPADECOEFFICIENTS returns coefficients of numerator
            %    of [M/M] Pade approximant, where M = 3,5,7,9,13.
            switch m
                case 3
                    c = [120, 60, 12, 1];
                case 5
                    c = [30240, 15120, 3360, 420, 30, 1];
                case 7
                    c = [17297280, 8648640, 1995840, 277200, 25200, 1512, 56, 1];
                case 9
                    c = [17643225600, 8821612800, 2075673600, 302702400, 30270240, ...
                         2162160, 110880, 3960, 90, 1];
                case 13
                    c = [64764752532480000, 32382376266240000, 7771770303897600, ...
                         1187353796428800,  129060195264000,   10559470521600, ...
                         670442572800,      33522128640,       1323241920,...
                         40840800,          960960,            16380,  182,  1];
            end
        end
    end
    
%%%%Nested Functions%%%%
end
