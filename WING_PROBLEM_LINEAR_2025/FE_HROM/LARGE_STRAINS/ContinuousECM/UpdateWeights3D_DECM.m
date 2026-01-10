function UpdateWeights3D_DECM

if nargin == 0
    load('tmp1.mat')
end

fhandle = guihandles(gcf) ; % --> Tag
hread = getfield(fhandle,'slide_tag');
VAL = get(hread,'Value') ;
iter = round(VAL)      ;




DATA = guidata(gcf);


 UpdateWeights3Dloc_DECM(DATA,iter) ; 


 