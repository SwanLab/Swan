function [Umixed_n,Smixed_n,Vmixed_n] =  MIXEDsvdVEL_DISP(Udisp,Sdisp,Vdisp,Uvel,Svel,Vvel,TIME_STEPS)


SVdisp_n = bsxfun(@times,Vdisp(TIME_STEPS,:)',Sdisp)  ;
SVvel_n = bsxfun(@times,Vvel(TIME_STEPS,:)',Svel)  ;
[UU,Smixed_n,Vmixed_n] = RSVDT([SVdisp_n;SVvel_n]) ;
Umixed_n = zeros(2*size(Udisp,1),length(Smixed_n));
ncolsDISP = 1:size(Udisp,2) ;
ncolsVEL = size(Udisp,2)+1:length(Smixed_n) ;
Umixed_n(1:size(Udisp,1),:) = Udisp*UU(ncolsDISP,:) ;
Umixed_n(size(Udisp,1)+1:end,:) = Uvel*UU(ncolsVEL,:) ;