function PlotPathPointsGauss1D(ELEMENTS_TO_PLOT,VAR_SMOOTH_FE,PATHiniPOINTS,xINITIAL,xOLD)

if nargin == 0
    load('tmp1.mat')
end


figure(15)
hold on
xlabel('x')
ylabel('y')

% Plot Elements affected during the iterations
% ---------------------------------------------
%         INDused = 1:length(POLYINFO.COEFFSpolynomial) ;
%         INDEMPT = cellfun(@isempty,POLYINFO.COEFFSpolynomial) ;
%         INDused(INDEMPT) = [];


% 
% INDused= unique(ELEMENTS_TO_PLOT) ;
% VAR_SMOOTH_FE  =DefaultField(VAR_SMOOTH_FE,'PlotContainingElementsFigureDECM_CECM',0) ;
% if VAR_SMOOTH_FE.PlotContainingElementsFigureDECM_CECM == 1
%     for ielem = 1:length(INDused)
%         elemLOC = INDused(ielem) ;
%         NODESloc  = VAR_SMOOTH_FE.CN(elemLOC,VAR_SMOOTH_FE.IND_POLYG_ELEMENT)' ;
%         COORloc = VAR_SMOOTH_FE.COOR(NODESloc,:) ;
%         plot(COORloc(:,1),COORloc(:,2),'LineWidth',1) ;
%     end
% end

% Plotting path during removal iterations
% PATHiniPOINTS

% for ipoint =1:size(PATHiniPOINTS,1)
%     locPATH = PATHiniPOINTS(ipoint,:) ;
%     indN = cellfun(@isempty,locPATH) ;
%     locPATH = locPATH(indN==0) ;
%     locPATH = cellfun(@transpose,locPATH,'UniformOutput',false) ;
%     COOR = cell2mat(locPATH) ;
%     %           vxlabels = arrayfun(@(n) {sprintf('k%d', n)}, (1:size(COOR,2))');
%     %Hpl = text(COOR(1,:), COOR(2,:), vxlabels,'HorizontalAlignment','center');
%     
%     plot(COOR(1,:),COOR(2,:),'Marker','.') ;
%     
%     
% end



%  end
h = plot(xINITIAL(:,1),xINITIAL(:,2),'o','Color',[0 0 1],'MarkerSize',8);
h2 = plot(xOLD(:,1),xOLD(:,2),'x','Color',[1 0 0],'MarkerSize',8);
LLL{1} = ['Before opt., m =',num2str(size(xINITIAL,1))] ;
LLL{2} = ['After opt., m =',num2str(size(xOLD,1))] ;
legend([h,h2],LLL);
axis equal