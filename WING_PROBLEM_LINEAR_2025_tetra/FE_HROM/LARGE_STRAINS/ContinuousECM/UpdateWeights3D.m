function UpdateWeights3D


fhandle = guihandles(gcf) ; % --> Tag
hread = getfield(fhandle,'slide_tag');
VAL = get(hread,'Value') ;
iter = round(VAL)      ;




DATA = guidata(gcf);


 UpdateWeights3Dloc(DATA,iter) ; 





% iter_before = max(1,iter-1) ;
% OldWeight = DATA.wALL(:,iter_before) ;
% %
% set(DATA.h_old,'Ydata',OldWeight) ;


% DATA.stepplot = round(Value);
% [DATA]= May29_step(DATA);
% guidata(gcf,DATA) ;
