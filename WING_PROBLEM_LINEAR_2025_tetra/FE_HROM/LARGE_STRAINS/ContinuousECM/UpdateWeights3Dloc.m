function  UpdateWeights3Dloc(DATA,iter)


NewWeight = DATA.wALL(:,iter) ;
NewPoint = DATA.xALL{iter} ; %(:,iter) ;

for iplot = 1:length(DATA.h)
    set(DATA.h(iplot),'Xdata',[NewPoint(iplot,1),NewPoint(iplot,1)]) ;
    set(DATA.h(iplot),'Ydata',[NewPoint(iplot,2),NewPoint(iplot,2)]) ;
    set(DATA.h(iplot),'Zdata',[NewPoint(iplot,3),NewPoint(iplot,3)]) ;
    
    wPOINT = NewWeight(iplot) ;
    if wPOINT == 0
        set(DATA.h(iplot),'Marker','none') ;
    else
        MarkerSizeLoc =  DATA.MarkerSizeMin +(DATA.MarkerSizeMax-DATA.MarkerSizeMin)/(DATA.wMAX-DATA.wMIN)*(wPOINT-DATA.wMIN) ;
        set(DATA.h(iplot),'MarkerSize',abs(MarkerSizeLoc)) ;
        
    end
    
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
    
    
    set(DATA.htitle,'String',['STEP =',num2str(iter),';  Number of points = ',num2str(npoints),' (of ',num2str(size(DATA.wALL,1)),')']) ;
    %  htitle = title(['First stage (elimination): STEP =',num2str(iter),';  Number of points = ',num2str(npoints),' (of ',num2str(size(wALL,1)),')'])  ;
    
    
    if iter == length(DATA.xALL)
        set(DATA.htitle,'String',['Final CECM rule with npoints = ',num2str(npoints),' (of ',num2str(size(DATA.wALL,1)),')']) ;
        
        delete(DATA.hPOINT)
        delete(DATA.legend) ; %,'String','Final CECM rule') ;%,'String','progressive reduction of weights)')
        
    else
        if iterNEXT > DATA.INDEXES_FIRST_STAGE(end)
            %  set(DATA.legend,'String','') ;%,
            % set(DATA.hPOINT,'Marker','none') ;
        else
            if iterNEXT <=  (DATA.niter)
                IndNEXT = find(DATA.wALL(:,iterNEXT) == 0) ;
                IndRemoveLOC = setdiff(IndNEXT,IndNOW) ;
                xRED = DATA.xALL{iter}(IndRemoveLOC,:) ;
                % wRED = DATA.wALL(IndRemoveLOC,iter) ;
                set(DATA.hPOINT,'Xdata',xRED(1)) ;
                set(DATA.hPOINT,'Ydata',xRED(2)) ;
                set(DATA.hPOINT,'Zdata',xRED(3)) ;
            else
                %set(DATA.hPOINT,'MarkerSize',0.1,'Marker','.') ;
            end
        end
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
    % istep = max(istep,2) ;
    ControlPointIndex =  DATA.CONTROL_POINTS(istep-1) ;
    
    
    set(DATA.htitle,'String',['2nd stage (progressive reduction); step=',num2str(iter),';Number of points =',num2str(npoints),' (of ',num2str(size(DATA.wALL,1)),')']) ;
    
    if iter == length(DATA.xALL)
        set(DATA.htitle,'String',['Final CECM rule with npoints = ',num2str(npoints),' (of ',num2str(size(DATA.wALL,1)),')']) ;
        
        delete(DATA.hPOINT)
        delete(DATA.legend) ; %,'String','Final CECM rule') ;%,'String','progressive reduction of weights)')
        
    else
        set(DATA.legend,'String','Controlled point') ;%,'String','progressive reduction of weights)')
        
        xRED = DATA.xALL{iter}(ControlPointIndex,:) ;
        %  wRED = DATA.wALL(ControlPointIndex,iter) ;
        set(DATA.hPOINT,'Xdata',xRED(1)) ;
        set(DATA.hPOINT,'Ydata',xRED(2)) ;
        set(DATA.hPOINT,'Zdata',xRED(3)) ;
    end
    
    
end
