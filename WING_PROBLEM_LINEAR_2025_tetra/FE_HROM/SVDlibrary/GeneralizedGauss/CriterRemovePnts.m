function Scrit = CriterRemovePnts(DATALOC,wNEW,PHIk_y,b,PHI_der,xNEW,xINT)

%CRITER_REMOVE = DATALOC.CRITERION_ELIMINATION;
%if CRITER_REMOVE==0
    Scrit = wNEW.*sum(PHIk_y.*PHIk_y,2) ;
%elseif   CRITER_REMOVE==1
%    Scrit = wNEW ;
%elseif  CRITER_REMOVE==2
 %   error('This option is not efficient... Set = 0')
    % Full Jacobian
  %  bNEW = PHIk_y'*wNEW ; 
%     F = b-bNEW ;
%     D= JAcobianGaussInt2D(PHI_der,xNEW,wNEW,PHIk_y,xINT) ;
%     %
%     indALL = 1:size(D,2) ; 
%     Scrit = zeros(size(wNEW));
%     for  ipoint = 1:length(wNEW)
%         indP = setdiff(indALL,ipoint:length(wNEW):size(D,2)) ;
%          
%         dq = D(:,indP)'*((D(:,indP)*D(:,indP)')\F) ; 
%         Scrit(ipoint) = norm(dq) ; 
%     end
% end