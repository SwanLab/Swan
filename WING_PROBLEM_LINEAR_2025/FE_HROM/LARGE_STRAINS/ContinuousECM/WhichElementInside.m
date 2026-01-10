function elemCONTAINER = WhichElementInside(xLOC,INDnear,VAR_SMOOTH_FE,IND_POLYG) ; 

ielem = 1;
elemCONTAINER = [] ;
   [ELEMnear, aaaa ]= find(VAR_SMOOTH_FE.CN == INDnear) ;  % Elements sharing INDnear.

while ielem <= length(ELEMnear)
    elemLOC = ELEMnear(ielem) ;
    
    [inELEM,onELEM] = IsInsideALL(xLOC,VAR_SMOOTH_FE,elemLOC,IND_POLYG)  ;
    
    if inELEM == 1 || onELEM == 1
        elemCONTAINER  = elemLOC ;
        %             if elemCONTAINER == 99
        %                 disp('Borrar esto')
        %             end
        break
    end
    ielem = ielem + 1;
end