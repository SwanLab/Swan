function  UpdateWeights3Dloc_DECM(DATA,iter)

if nargin == 0
    load('tmp1.mat')
end

NewWeight = DATA.wALL(:,iter) ;
%NewPoint = DATA.xALL{iter} ; %(:,iter) ;

for iplot = 1:length(DATA.h)
   
    
    wPOINT = NewWeight(iplot) ;
    if wPOINT == 0
        set(DATA.h(iplot),'Marker','none') ;
    else
        MarkerSizeLoc =  DATA.MarkerSizeMin +(DATA.MarkerSizeMax-DATA.MarkerSizeMin)/(DATA.wMAX-DATA.wMIN)*(wPOINT-DATA.wMIN) ;
        set(DATA.h(iplot),'Marker','.','MarkerSize',abs(MarkerSizeLoc)) ;
        
    end
    
end
 


%iterNEXT = iter + 1;
IndNOW = find(DATA.wALL(:,iter) == 0) ;
npoints = size(DATA.wALL,1)-length(IndNOW) ;


%htitle = title(['DECM: error (%) =',num2str(errorDECM(iter)),';  Number of points = ',num2str(npoints),' (of ',num2str(size(wALL,1)),')'])  ;


STRING_txt = ['DECM:  error (%)  =',num2str(DATA.errorDECM(iter)*100),';  Number of points = ',num2str(npoints),' (of ',num2str(size(DATA.wALL,1)),')']  ;

set(DATA.htitle,'String',STRING_txt) ;
 