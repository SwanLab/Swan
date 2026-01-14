function [xOLD,wOLD,DATALOC,POINTS_all,iter] = ...
    Iter_RemPnt_AtTheEnd(xOLD,wOLD,xINT,PHI,PHI_der,b,DATALOC) ;

if nargin == 0 
    load('tmp.mat')
end
CONVERGENCE = 1 ;
iter = 1;

while CONVERGENCE ==1 & length(xOLD)>1
    %Iter_RemPnt2Dapprox_END --> Copy of Iter_RemPnt2Dapprox
    [xNEW,wNEW,CONVERGENCE,DATALOC] = Iter_RemPnt2Dapprox_END(xOLD,wOLD,xINT,PHI,PHI_der,b,DATALOC) ;
    if CONVERGENCE == 1
        xOLD = xNEW ;
        wOLD = wNEW ;
        POINTS_all.x{iter} = xNEW ;
        POINTS_all.w{iter} = wNEW ;
        iter = iter + 1 ;
    end
end