function   INDICES = DEIMfun(PHI,NPOINTS)
%---------------------------------------------------
% See  --> Greedy_Hierarchical_cycle.m  ------------
%---------------------------------------------------
if nargin == 0
    
    load('tmp.mat')
end
if nargin == 1
    NPOINTS = size(PHI,2) ;
end

if NPOINTS < size(PHI,2)
    error('Inconsistent number of points')
end

Uk = PHI(:,1) ;
m = size(PHI,2);
INDICES = zeros(m,1) ;
k = 1 ;

[rho ind] = max(abs(PHI(:,k)));
INDICES(k) = ind ;

%dbstop('156')
%  hhhh = waitbar(0,'Computing norm(max(abs(r))) ...');

for k = 2:m
    %     waitbar(k/m);
    ind_k = INDICES(1:(k-1));
    %U = PHI(:,1:(k-1));
    U_red = Uk(ind_k,:);
    c = U_red\PHI(ind_k,k);
    r = PHI(:,k)-Uk*c ;
    [  rho_max   ind_kp1   ] = max(abs(r));
    Uk = [Uk PHI(:,k)];
    INDICES(k) = ind_kp1 ;
end

npointTOT= size(PHI,1);
RAND_POINTS = randi(npointTOT,1,npointTOT);

if NPOINTS >m
    if NPOINTS > npointTOT
        INDICES = 1: npointTOT ;
    else
    % Selecting additional entries. But how ?
    %
    
    ENCONTRADOS = zeros(npointTOT,1) ;
    ENCONTRADOS(INDICES) = 1 ;
    nrest = NPOINTS-m ;
    condMIN = 50000000 ;
  
    for irest = 1:nrest
        ipoint = m + irest ;
        inewLOC = 1
        
        while inewLOC <=npointTOT
            inew = RAND_POINTS(inewLOC) ;
            if ENCONTRADOS(inew)==0
                INDICES_NEW = [INDICES; inew];
                PHI_RED = PHI(INDICES_NEW,:);
                M = PHI_RED'*PHI_RED ;
                if 1/rcond(M)<condMIN
                    INDICES =  INDICES_NEW  ;
                    ENCONTRADOS(inew) = 1 ;
                    break
                end
            else
                inewLOC = inewLOC+1 ;
            end
        end
    end
    end
end

