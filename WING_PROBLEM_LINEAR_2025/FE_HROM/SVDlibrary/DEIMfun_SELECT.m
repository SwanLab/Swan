function   IND_SNAP = DEIMfun_SELECT(SNAP,DATAoffline)
%---------------------------------------------------
%---------------------------------------------------
if nargin == 0
    load('tmp1.mat')
end
 DATAoffline = DefaultField(DATAoffline,'TOLERANCE_SELECTION_DEIM',1e-6) ; 
  DATAoffline = DefaultField(DATAoffline,'NCLUSTERS_BASIS_DISP',[]) ; 

% First snapshot 
% --------------
TOL_convergence = 1e-12; 

nsnap = size(SNAP,2) ; 

SNAP_n = sqrt(sum(SNAP.*SNAP,1)) ;

% We remove from the candidate set those snapshots with norm below TOL_convergence

TOL_REMOVE = 0.01*TOL_convergence ; 
IND_CANDIDATES=  find( SNAP_n >TOL_REMOVE) ; 

% 

SNAP = bsxfun(@times,SNAP',1./SNAP_n')' ; % Normalization 

[~,IND_SNAP] =max(SNAP_n) ; 
IND_SNAP = IND_SNAP(1) ; 

iLOC=1 ; 
rNmin = 1e20 ; 

if isempty(DATAoffline.NCLUSTERS_BASIS_DISP)
    ncandidates= length(IND_CANDIDATES) ; 
else
    ncandidates = DATAoffline.NCLUSTERS_BASIS_DISP-1 ; 
    TOL_convergence = 0; 
end

while iLOC<=ncandidates && rNmin>TOL_convergence
    
    % Projection of SNAP onto the span of SPAN(:,IND_SNAP)
    c = SNAP(:,IND_SNAP)\SNAP(:,IND_CANDIDATES) ;
    % We select as next index the one with lower c, as this indicates
    r = SNAP(:,IND_CANDIDATES)- SNAP(:,IND_SNAP)*c ;
    rN = sqrt(sum(r.*r,1)) ;
    [rNmin,IIIminLOC] = max(rN) ;
    IIImin = IND_CANDIDATES(IIIminLOC) ; 
    disp(['Iter = ',num2str(iLOC),'; SNAP = ',num2str(IIImin),'; RES =',num2str(rNmin)]) ; 
    
    
    IND_SNAP = [IND_SNAP,IIImin] ;
    
    iLOC = iLOC +1 ;
    
    
end

disp(['SELECTED snapshots = ',num2str(sort(IND_SNAP))])



% Uk = PHI(:,1) ;
% m = size(PHI,2);
% INDICES = zeros(m,1) ;
% k = 1 ;
% 
% [rho ind] = max(abs(PHI(:,k)));
% INDICES(k) = ind ;
% 
% disp(['iter = ',num2str(k),'; SNAPSHOT =',num2str(ind)])
% 
% %dbstop('156')
% %  hhhh = waitbar(0,'Computing norm(max(abs(r))) ...');
% 
% for k = 2:m
%     %     waitbar(k/m);
%     ind_k = INDICES(1:(k-1));
%     %U = PHI(:,1:(k-1));
%     U_red = Uk(ind_k,:);
%     c = U_red\PHI(ind_k,k);
%     r = PHI(:,k)-Uk*c ;
%     [  rho_max   ind_kp1   ] = max(abs(r));
%     
%         disp(['iter = ',num2str(k),'; SNAPSHOT =',num2str(ind_kp1),' NORM RES =',num2str(norm(r))])
% 
%     
%     Uk = [Uk PHI(:,k)];
%     INDICES(k) = ind_kp1 ;
% end
% 
% npointTOT= size(PHI,1);
% RAND_POINTS = randi(npointTOT,1,npointTOT);
% 
% if NPOINTS >m
%     if NPOINTS > npointTOT
%         INDICES = 1: npointTOT ;
%     else
%     % Selecting additional entries. But how ?
%     %
%     
%     ENCONTRADOS = zeros(npointTOT,1) ;
%     ENCONTRADOS(INDICES) = 1 ;
%     nrest = NPOINTS-m ;
%     condMIN = 50000000 ;
%   
%     for irest = 1:nrest
%         ipoint = m + irest ;
%         inewLOC = 1
%         
%         while inewLOC <=npointTOT
%             inew = RAND_POINTS(inewLOC) ;
%             if ENCONTRADOS(inew)==0
%                 INDICES_NEW = [INDICES; inew];
%                 PHI_RED = PHI(INDICES_NEW,:);
%                 M = PHI_RED'*PHI_RED ;
%                 if 1/rcond(M)<condMIN
%                     INDICES =  INDICES_NEW  ;
%                     ENCONTRADOS(inew) = 1 ;
%                     break
%                 end
%             else
%                 inewLOC = inewLOC+1 ;
%             end
%         end
%     end
%     end
% end
% 
