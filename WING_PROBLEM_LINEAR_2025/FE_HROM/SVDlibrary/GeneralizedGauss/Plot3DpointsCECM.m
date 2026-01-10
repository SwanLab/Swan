function Plot3DpointsCECM(DATALOC,xDECM,xCECM,xINITIAL,wCECM,wDECM,zCECM,DATAadd)

if nargin ==0
    load('tmp.mat')
end

figure(15)
hold on
xlabel('x')
ylabel('y')
zlabel('z')
grid on
%     if DATALOC.INITIAL_DISCRETIZATION_TENSORPRODUCT ==1
%         plot3(xINITIAL(:,1),xINITIAL(:,2),xINITIAL(:,3),'.','Color',[0 1 0])
%     end
h = plot3(xDECM(:,1),xDECM(:,2),xDECM(:,3),'o','Color',[0 0 1],'MarkerSize',8);
h2 = plot3(xCECM(:,1),xCECM(:,2),xCECM(:,3),'x','Color',[1 0 0],'MarkerSize',8);
hseed = plot3(xINITIAL(zCECM,1),xINITIAL(zCECM,2),xINITIAL(zCECM,3),'*','Color',[1 0 1],'MarkerSize',8);
LLL{1} = ['Before opt., m =',num2str(size(xDECM,1))] ;
LLL{2} = ['After opt., m =',num2str(size(xCECM,1))] ;
LLL{3} = ['Seed points '] ;
legend([h,h2,hseed],LLL);




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

     figure(16)
    hold on
    bar(sort(wCECM))
    ylabel('Weights (after optimization)')
    figure(17)
    hold on
    bar(sort(wDECM))
    ylabel('Weights (before optimization)')
    
 
%diary off
