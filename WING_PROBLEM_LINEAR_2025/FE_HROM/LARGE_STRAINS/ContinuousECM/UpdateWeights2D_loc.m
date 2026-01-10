function UpdateWeights2D


fhandle = guihandles(gcf) ; % --> Tag
hread = getfield(fhandle,'slide_tag');
VAL = get(hread,'Value') ;
iter = round(VAL)      ;




DATA = guidata(gcf);
NewWeight = DATA.wALL(:,iter) ;
NewPoint = DATA.xALL{iter} ; %(:,iter) ;

for iplot = 1:length(DATA.h)
    set(DATA.h(iplot),'Zdata',[0,NewWeight(iplot)]) ;
    set(DATA.h(iplot),'Xdata',[NewPoint(iplot,1),NewPoint(iplot,1)]) ;
    set(DATA.h(iplot),'Ydata',[NewPoint(iplot,2),NewPoint(iplot,2)]) ;
end

IS_FIRST_STAGE = 1;
if ~isempty( DATA.INDEXES_FIRST_STAGE)
    if  iter > DATA.INDEXES_FIRST_STAGE(end)
        IS_FIRST_STAGE = 0 ;
    end
end

if DATA.Include2ndStageIterations_PlotEvolutionWeights == 0 || IS_FIRST_STAGE
    
    
    
    
    iterNEXT = iter + 1;
    IndNOW = find(DATA.wALL(:,iter) == 0) ;
    npoints = size(DATA.wALL,1)-length(IndNOW) ;
    
    
    set(DATA.htitle,'String',['Number of points = ',num2str(npoints),' (of ',num2str(size(DATA.wALL,1)),')']) ;
    
    if iterNEXT <=  (DATA.niter)
        IndNEXT = find(DATA.wALL(:,iterNEXT) == 0) ;
        IndRemoveLOC = setdiff(IndNEXT,IndNOW) ;
        xRED = DATA.xALL{iter}(IndRemoveLOC,:) ;
        % wRED = DATA.wALL(IndRemoveLOC,iter) ;
        set(DATA.hPOINT,'Xdata',xRED(1)) ;
        set(DATA.hPOINT,'Ydata',xRED(2)) ;
    else
        %set(DATA.hPOINT,'MarkerSize',0.1,'Marker','.') ;
    end
    
else
    
    
    
    
    %   iterNEXT = iter + 1;
    IndNOW = find(DATA.wALL(:,iter) == 0) ;
    npoints = size(DATA.wALL,1)-length(IndNOW) ;
    
    %  ControlPointIndex  = DATA.CONTROL_POINTS;
    % IndexesSecondStage = DATA.INDEXES_SECOND_STAGE;
    istep  = 1;   ISITER = 0 ;
    while istep <= length(DATA.INDEXES_SECOND_STAGE) &&  all(ISITER==0)
        ISITER=  ismember(DATA.INDEXES_SECOND_STAGE{istep},iter) ;
        istep = istep +1;
    end
    ControlPointIndex =  DATA.CONTROL_POINTS(istep-1) ;
    
    set(DATA.legend,'String','Controlled point') ;%,'String','progressive reduction of weights)')
    
    %     if ~isempty(DATA.LEGEND_GAUSS)
    %     DATA.legend = legend([hPOINT,DATA.h_cecm,DATA.h_gauss],{'Point to be eliminated','CECM rule',DATA.LEGEND_GAUSS})
    % else
    % DATA.legend = legend([hPOINT,DATA.h_cecm],{'Point to be eliminated','CECM rule'})
    % end
    %
    
    set(DATA.htitle,'String',['2nd stage (progressive reduction): Number of points = ',num2str(npoints),' (of ',num2str(size(DATA.wALL,1)),')']) ;
    
    %     if iterNEXT <=  (DATA.niter)
    % IndNEXT = find(DATA.wALL(:,iterNEXT) == 0) ;
    %         IndRemoveLOC = setdiff(IndNEXT,IndNOW) ;
    xRED = DATA.xALL{iter}(ControlPointIndex,:) ;
    %  wRED = DATA.wALL(ControlPointIndex,iter) ;
    set(DATA.hPOINT,'Xdata',xRED(1)) ;
    set(DATA.hPOINT,'Ydata',xRED(2)) ;
    %     else
    %        set(DATA.hPOINT,'MarkerSize',0.1,'Marker','.') ;
    % end
    
end





% iter_before = max(1,iter-1) ;
% OldWeight = DATA.wALL(:,iter_before) ;
% %
% set(DATA.h_old,'Ydata',OldWeight) ;


% DATA.stepplot = round(Value);
% [DATA]= May29_step(DATA);
% guidata(gcf,DATA) ;
