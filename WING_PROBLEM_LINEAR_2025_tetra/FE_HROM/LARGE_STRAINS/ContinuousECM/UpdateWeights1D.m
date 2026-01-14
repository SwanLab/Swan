function UpdateWeights1D

DATA = guidata(gcf);
fhandle = guihandles(gcf) ; % --> Tag
hread = getfield(fhandle,'slide_tag');
VAL = get(hread,'Value') ;
iter = round(VAL)      ;

UpdateWeights1Dloc(DATA,iter) ; 