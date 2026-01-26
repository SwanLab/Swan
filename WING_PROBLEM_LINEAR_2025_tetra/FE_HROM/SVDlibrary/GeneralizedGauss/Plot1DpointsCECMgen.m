function Plot1DpointsCECMgen(DATALOC,xDECM,xCECM,xINITIAL,wCECM,wDECM,zCECM,DATAadd,DATAapprox,zDECM)

if nargin ==0
    load('tmp.mat')
end

figure(15)
hold on
xlabel('x')
ylabel('weight')

grid on
%     if DATALOC.INITIAL_DISCRETIZATION_TENSORPRODUCT ==1
%         plot3(xINITIAL(:,1),xINITIAL(:,2),xINITIAL(:,3),'.','Color',[0 1 0])
%     end
h(1) = bar(xDECM,wDECM);
h(2) = bar(xCECM,wCECM);
%h(3) = plot(xINITIAL(zCECM,1),xINITIAL(zCECM,2),'*','Color',[1 0 1],'MarkerSize',8);
LLL{1} = ['DECM, m =',num2str(size(xDECM,1)),' (error =',num2str(DATAadd.errorDECM),' %)'] ;
LLL{2} = ['CECM, m =',num2str(size(xCECM,1)),' (error =',num2str(DATAadd.errorCECM),' %)'] ;
%LLL{3} = ['Seed points '] ;


if isfield(DATAadd,'wGAUSS')
    h(3)  = bar(DATAadd.xGAUSS,DATAadd.wGAUSS) ;
    LLL{3} = ['Gaussian rule  m =',num2str(size(DATAadd.xGAUSS,1)),' points (error = ',num2str(DATAadd.errorGAUSS),', %)' ] ;
end
% 
% DATALOC = DefaultField(DATALOC,'PLOT_HISTORY_OPTIMAL_POINTS',0) ; 
%  
% if  DATALOC.PLOT_HISTORY_OPTIMAL_POINTS == 1
%       zHISTORY = DATAapprox.INFO_iterations.z ; 
%       xHISTORY = DATAapprox.INFO_iterations.x ; 
%       xHISTORY_zCECM = zeros(length(zCECM),length(zHISTORY)) ; 
%       yHISTORY_zCECM = xHISTORY_zCECM ; 
%     for kiter =1:length(zHISTORY)
%         [AA,IND_int] = intersect(zHISTORY{kiter},zCECM) ; 
%         xLOC = xHISTORY{kiter}(IND_int,:) ; 
%         xHISTORY_zCECM(:,kiter) = xLOC(:,1)' ; 
%         yHISTORY_zCECM(:,kiter) = xLOC(:,2)' ; 
%         
%     end
%     
%     hold on 
%     for iii = 1:length(zCECM) 
%     h(end+1) = plot(xHISTORY_zCECM(iii,:),yHISTORY_zCECM(iii,:),'Marker','.') ; 
%     LLL{end+1} = ['Traj. point = ',num2str(zCECM(iii))] ; 
%     end
%     
%     
%     
% end


legend(h,LLL);




%
%
% % Nearest points
% % ----------------
% STUDY_NEAREST_POINT = 1 ;
% if STUDY_NEAREST_POINT ==1
%     dt = DelaunayTri(xINITIAL(:,1),xINITIAL(:,2),xINITIAL(:,3));
%     NP = nearestNeighbor(dt, xCECM(:,1),xCECM(:,2),xCECM(:,3));
%
%     xNEAR = xINITIAL(NP,:); % Nearest point
%     EVALUATE_GRADIENT = 0 ;
%     VSinv = [] ;
%     g = Evaluate_Basis_Grad_Analytic(xNEAR,VSinv,DATALOC,EVALUATE_GRADIENT)  ;
%
%
%
%     % Weights
%     wNEAR = g\b ;
%     bNEAR = PHIk_y'*wNEAR ;
%     errorINTnear = norm(bNEAR-b)/norm(b)*100;
%     disp(['Error nearest points (%) =',num2str(errorINT)])
%     %%%
% end
% 
% SHOW_AFTER_BEFORE_SORTED = 0 ; 
% 
% if SHOW_AFTER_BEFORE_SORTED == 0
%     figure(16) 
%     subplot(2,1,1) 
%     hold on 
%     xlabel('Index   point (DECM, local)')
%     ylabel('DECM weights')
%     bar(wDECM)
%     
%     % Iteration at which each point is eliminated
%     % --------------------------------------------
%     DATALOC = DefaultField(DATALOC,'SHOW_ITERATION_WEIGHTS_ARE_ELIMINATED',0) ; 
%     if DATALOC.SHOW_ITERATION_WEIGHTS_ARE_ELIMINATED == 1
%     zHISTORY = DATAapprox.INFO_iterations.z ; 
%     zREMAIN= zDECM ; 
%     for kiter = 1:length(zHISTORY)
%      [ z_elim ]=   setdiff(zREMAIN,zHISTORY{kiter}) ;
%        AA = find(zDECM == z_elim);  
%      yTEXT = wDECM(AA)+0.1 ; 
%      text(AA,yTEXT,['k=',num2str(kiter)]) ; 
%      zREMAIN = zHISTORY{kiter} ; 
%     end
%     end
%     
%     
%     
%     % -------------------------------------------
%     
%     subplot(2,1,2)
%     hold on 
%      xlabel('Index   point (DECM, local)')
%     ylabel('CECM weights')
%     
%     [II,JJ,KK] = intersect(zDECM,zCECM) ;
%     wCECMsparse = zeros(size(wDECM)) ; 
%     wCECMsparse(JJ) = wCECM ; 
%     
%     
%       bar(wCECMsparse)
% else
%     
%     figure(16)
%     hold on
%     bar(sort(wCECM))
%     ylabel('Weights (after optimization)')
%     figure(17)
%     hold on
%     bar(sort(wDECM))
%     ylabel('Weights (before optimization)')
%     
% end
% 
% 
%  
% 
% %diary off
