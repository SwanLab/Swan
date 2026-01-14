function [xNEW,wNEW,CONVERGENCE,DATALOC,POLYINFO] = Iter_RemPnt2Dapprox(xNEW,wNEW,xINT,PHI,PHI_der,b,DATALOC,VAR_SMOOTH_FE)

%dbstop('5')
if nargin == 0
    load('tmp1.mat')
    
    %   load('tmpPR.mat')
    %    xNEW = x_ProducRule
    %    wNEW = w_ProducRule
    
    load('tmpMINE.mat','wGAUSS','xGAUSS')
    xNEW = xGAUSS ;
    wNEW = wGAUSS ;
    
end

DATALOC = DefaultField(DATALOC,'MSGPRINT',{}) ;
% Evaluation of the basis functions at the points xNEW
 DATALOC.REMOVED_INDEX = [] ; 

% if isempty(DATALOC.FUN_GRADIENT_DIRECTLY)
%     PHIk_y = zeros(size(xNEW,1),size(PHI,2));
%
%     for i=1:size(PHI,2)
%         PHIloc = reshape(PHI(:,i),size(xINT{1},1),[]);
%         PHIk_y(:,i) = interp2(xINT{1},xINT{2},PHIloc,xNEW(:,1),xNEW(:,2),'cubic') ;
%     end
% else
% dbstop('20')
%     if DATALOC.APPROX_FUN__DERI.ACTIVE == 0
%         PHIk_y = EvaluateBasisFunctionAtX(xNEW,VSinv,DATALOC,0,1)  ;
%     else
% Approximated method, fun. and derivatives
% ------------------------------------------- April-2020
% disp('Approx. fun. and der. ...')
%  tic
%  PHIk_y_old = EvaluateBasisFunctionAtX(xNEW,VSinv,DATALOC,0,1)  ;

POLYINFO = [] ; % Information about local searchs (for FE_interpolation option)
[PHIk_y,~,POLYINFO] = EvaluateBasisFunctionAtX_approx(xNEW, DATALOC.APPROX_FUN__DERI,VAR_SMOOTH_FE,POLYINFO)  ;


% for ipoint = 1:length(DATALOC.PATHiniPOINTS) 
%     DATALOC.PATHiniPOINTS{ipoint} = xNEW(ipoint,:) ;  % Coordinates of the initial points at each iteration 
% end
% DATALOC.IND_POINTS_NOT_REMOVED = 1:length(DATALOC.PATHiniPOINTS)  ; 

%  EEE = norm(PHIk_y-PHIk_y_old,'fro')./norm(PHIk_y_old,'fro')*100 ;
%  disp(['Error approx = ',num2str(EEE),' %'])
%  toc
disp('...Done')
%   end
%end

errCOMP = norm(PHIk_y'*wNEW-b)./norm(b)     % This error measures how approximate is the SVD spatial decomp.

if DATALOC.iter == 1
    DATALOC.errorFITTING = errCOMP ;
    DATALOC.MSGPRINT{end+1}= '--------------------------------------------' ;
    DATALOC.MSGPRINT{end+1}= ['Integration appr. error at k= 0 -->',num2str(errCOMP*100),' %', '(It measures the quality in evaluating the function via   fitting)'];
end

% ------------------T
% Criterion for selecting the point to be removed: significance index
% This process is repeated
disp('*********************************************************************************************')
disp([' CUB. RULE  with    ', num2str(length(xNEW)-1),' points ( of ',num2str(DATALOC.mINI),')'])
disp('**********************************************************************************************')
TOL = DATALOC.TOL ;
% Criterion for removing points
S_crit = CriterRemovePnts(DATALOC,wNEW,PHIk_y,b,PHI_der,xNEW,xINT);

[MMM indREM ]= sort(S_crit) ;

xOLD = xNEW ;
wOLD = wNEW ;
% ----------------------------------------------------
% New set of points (and weights)
m = size(xNEW,1) ;
p = size(PHI,2) ;
% Generalized variable
iremove = 1 ;
SALIR = 0 ;
%dbstop('43')
normB = norm(b);
%dbstop('37')
while SALIR == 0 & iremove <=size(xOLD,1)
    %----------------------------------------------------
    % Remove point indREM(iremove)
    xNEW(indREM(iremove),:) = [] ;
    wNEW(indREM(iremove)) = [] ;     
    
    disp('----------------------------------------')
    disp(['Removing point =',num2str(indREM(iremove)),' (irem =',num2str(iremove),')'])
    disp('---------------------------------------')
    kMAX = DATALOC.kMAX  ;
    k=1;
    nF = 1e10 ;
    CONVERGENCE = 1;
    nF_old = nF ;
    nITER_increase = 0 ;
    maxITER_increase = 4 ;
    xORG = xNEW ;
    % wORG = wNEW ;    
    while  nF>=TOL & k<= kMAX & nITER_increase<=maxITER_increase       
        
        [xkp1, wkp1, nF, ISNEGATIVE,POLYINFO  ] = IterNR_GAuss2Dapprox(wNEW,PHI,b,xNEW,PHI_der,xINT,DATALOC,VAR_SMOOTH_FE,POLYINFO) ;
        nF = nF/normB ;
        
        DATALOC = DefaultField(DATALOC,'hx',[]) ;
        dX = abs(xkp1-xORG) ;
        if ~isempty(DATALOC.hx)
            dxREL = max(dX(:,1))/DATALOC.hx*100 ;
            dyREL = max(dX(:,2))/DATALOC.hy*100 ;
            maxDD = max(DATALOC.dxRELmax , DATALOC.dyRELmax )  ;
            disp(['Iteration k=',num2str(k),',  error =',num2str(nF),', dxMAX(%)=',num2str(dxREL), ', dyMAX(%)=',num2str(dyREL),...
                ' MAXdxy(%) =',num2str(maxDD)])
        else
            disp(['Iteration k=',num2str(k),',  error =',num2str(nF)])
        end
        
        if ISNEGATIVE ==1
            break
        end       
        if nF>nF_old
            nITER_increase = nITER_increase + 1 ;
        end
        %%% NEW ITERATIOn
        %  disp('Iteration k=2')
        xNEW = xkp1;
        wNEW = wkp1 ;        
        k = k+1 ;
        nF_old = nF ;
    end
    % ------------------------------------------------------
    
    if ISNEGATIVE ==1 | (k>kMAX  &  nF>DATALOC.TOL_low ) | nITER_increase >maxITER_increase        
        iremove = iremove + 1;
        xNEW = xOLD ;
        wNEW = wOLD ;        
    else
        SALIR = 1 ;         
    end
end

if  iremove > length(xOLD)
    CONVERGENCE = 0 ;
    
else
    DATALOC.REMOVED_INDEX = indREM(iremove) ; 
     
end
