function UpdateWeights2D


fhandle = guihandles(gcf) ; % --> Tag
hread = getfield(fhandle,'slide_tag');
VAL = get(hread,'Value') ;
iter = round(VAL)      ;




DATA = guidata(gcf);

UpdateWeights2Dloc(DATA,iter)