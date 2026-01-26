%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2012
%--------------------------------------------------------------------------

function Cp = thrC(C,ro)

if (nargin < 2)
    ro = 1;
end

if (ro < 1)
    N = size(C,2);
    Cp = zeros(N,N);
    [S,Ind] = sort(abs(C),1,'descend');
    for i = 1:N
        cL1 = sum(S(:,i));
        stop = false;
        cSum = 0; t = 0;
        while (~stop)
            t = t + 1;
            cSum = cSum + S(t,i);
            if ( cSum >= ro*cL1 )
                stop = true;
                Cp(Ind(1:t,i),i) = C(Ind(1:t,i),i);
            end
        end
    end
else
    Cp = C;
end