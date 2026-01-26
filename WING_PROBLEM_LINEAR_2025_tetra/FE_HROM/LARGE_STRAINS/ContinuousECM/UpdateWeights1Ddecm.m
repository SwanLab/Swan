function UpdateWeights1Ddecm

DATA = guidata(gcf);
fhandle = guihandles(gcf) ; % --> Tag
hread = getfield(fhandle,'slide_tag');
VAL = get(hread,'Value') ;
iter = round(VAL)      ;

UpdateWeights1DlocDECM(DATA,iter) ; 

% if iter == size(DATA.wALL,2)
%     
% end