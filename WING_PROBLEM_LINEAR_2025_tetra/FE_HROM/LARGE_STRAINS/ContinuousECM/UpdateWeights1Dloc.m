function UpdateWeights1Dloc(DATA,iter)


NewWeight = DATA.wALL(:,iter) ;
NewPoint = DATA.xALL(:,iter) ;

for iplot = 1:length(DATA.h)
    set(DATA.h(iplot),'Ydata',[0,NewWeight(iplot)]) ;
    set(DATA.h(iplot),'Xdata',[NewPoint(iplot),NewPoint(iplot)]) ;
end

IS_FIRST_STAGE = 1;
if iter > DATA.INDEXES_FIRST_STAGE(end)
    IS_FIRST_STAGE = 0 ;
end

if DATA.Include2ndStageIterations_PlotEvolutionWeights == 0 || IS_FIRST_STAGE
    
    
    
    
    iterNEXT = iter + 1;
    IndNOW = find(DATA.wALL(:,iter) == 0) ;
    npoints = size(DATA.wALL,1)-length(IndNOW) ;
    
    
    set(DATA.htitle,'String',['Number of points = ',num2str(npoints),' (of ',num2str(size(DATA.wALL,1)),')']) ;
    
    if iterNEXT <=  (DATA.niter)
        IndNEXT = find(DATA.wALL(:,iterNEXT) == 0) ;
        IndRemoveLOC = setdiff(IndNEXT,IndNOW) ;
        xRED = DATA.xALL(IndRemoveLOC,iter) ;
        wRED = DATA.wALL(IndRemoveLOC,iter) ;
        set(DATA.hPOINT,'Xdata',xRED) ;
    else
        %set(DATA.hPOINT,'MarkerSize',0.1,'Marker','.') ;
    end
    
else
    
    
    
    
    %   iterNEXT = iter + 1;
    IndNOW = find(DATA.wALL(:,iter) == 0) ;
    npoints = size(DATA.wALL,1)-length(IndNOW) ;
    
    %  ControlPointIndex  = DATA.CONTROL_POINTS;
    % IndexesSecondStage = DATA.INDEXES_SECOND_STAGE;
    istep  = 1; SALIR = 0 ; ISITER = 0 ;
    while istep <= length(DATA.INDEXES_SECOND_STAGE) &&  all(ISITER==0)
        ISITER=  ismember(DATA.INDEXES_SECOND_STAGE{istep},iter) ;
        istep = istep +1;
    end
    ControlPointIndex =  DATA.CONTROL_POINTS(istep-1) ;
    
    set(DATA.legend,'String','Controlled point') ;%,'String','progressive reduction of weights)')
    set(DATA.htitle,'String',['2nd stage (progressive reduction): Number of points = ',num2str(npoints),' (of ',num2str(size(DATA.wALL,1)),')']) ;
    
    %     if iterNEXT <=  (DATA.niter)
    % IndNEXT = find(DATA.wALL(:,iterNEXT) == 0) ;
    %         IndRemoveLOC = setdiff(IndNEXT,IndNOW) ;
    xRED = DATA.xALL(ControlPointIndex,iter) ;
    wRED = DATA.wALL(ControlPointIndex,iter) ;
    set(DATA.hPOINT,'Xdata',xRED) ;
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
