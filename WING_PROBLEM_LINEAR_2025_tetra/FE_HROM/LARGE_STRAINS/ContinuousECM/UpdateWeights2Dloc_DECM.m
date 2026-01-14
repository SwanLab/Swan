function  UpdateWeights2Dloc_DECM(DATA,iter)

if nargin == 0
    load('tmp1.mat')
end

NewWeight = DATA.wALL(:,iter) ;
%NewPoint = DATA.xALL{iter} ; %(:,iter) ;

for iplot = 1:length(DATA.h)
   
    
   % wPOINT = NewWeight(iplot) ;
    
    
     %h(iplot) = plot3(xALL(iplot,1)*ones(1,2) ,xALL(iplot,2)*ones(1,2),[0,wALL(iplot,1)],...
     %   'Color',[0,0,0],'LineWidth',AUXVAR.LineWidth);
    
%     if wPOINT == 0
%         set(DATA.h(iplot),'Marker','none') ;
%     else
%        MarkerSizeLoc =  DATA.MarkerSizeMin +(DATA.MarkerSizeMax-DATA.MarkerSizeMin)/(DATA.wMAX-DATA.wMIN)*(wPOINT-DATA.wMIN) ;
      %   set(DATA.h(iplot),'Marker','.','MarkerSize',abs(MarkerSizeLoc)) ;
      if  NewWeight(iplot) > 0
          set(DATA.h(iplot),'Zdata',[0,NewWeight(iplot)],'Marker','.','MarkerSize',6) ;
      else
            set(DATA.h(iplot),'Zdata',[0,NewWeight(iplot)]) ;
      end
%         
%     end
    
end
 


%iterNEXT = iter + 1;
IndNOW = find(DATA.wALL(:,iter) == 0) ;
npoints = size(DATA.wALL,1)-length(IndNOW) ;


%htitle = title(['DECM: error (%) =',num2str(errorDECM(iter)),';  Number of points = ',num2str(npoints),' (of ',num2str(size(wALL,1)),')'])  ;

DATA = DefaultField(DATA,'NumberOfPointsTotal',[]) ; 
if ~isempty(DATA.NumberOfPointsTotal)
STRING_txt = ['DECM:  error (%)  =',num2str(DATA.errorDECM(iter)*100),';  Number of points = ',num2str(npoints),' (of ',num2str((DATA.NumberOfPointsTotal)),')']  ;
else
    STRING_txt = ['DECM:  error (%)  =',num2str(DATA.errorDECM(iter)*100),';  Number of points = ',num2str(npoints)]  ;

end
set(DATA.htitle,'String',STRING_txt) ;
 