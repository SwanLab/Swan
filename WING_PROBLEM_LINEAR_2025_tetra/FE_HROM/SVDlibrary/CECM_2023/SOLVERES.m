function [xNEW,wNEW,VARCnew,SALIRloc,SALIR,DATALOC ] =...
    SOLVERES(normB,iremovLOC,wNEW,DATALOC,TOL,b,xNEW,VAR_SMOOTH_FE,...
    POLYINFO,VARCnew,iremove)

if nargin == 0
    load('tmp.mat')
end 
POINTS_F = [VARCnew.POINTSl(:); VARCnew.POINTSRp(:)];
ipoint_control = POINTS_F(iremovLOC) ;  % Index of the constrained point (global)
% --------------------------------------------------------------------------------------
% We would have to guess now whether this point belongs to VARCnew.POINTSl or
% VARCnew.POINTSRp. However, in this strategy, VARCnew.POINTSRp is always
% empty. Therefore, we simply make
%   if iremovLOC <= length(VARCnew.POINTSl)
VARCnew.POINTSl(iremovLOC)=[]  ;  % Remove the point from the unconstrained set
%POINTS_F=[VARCnew.POINTSl(:)'; VARCnew.POINTSRp(:)' ];  ;  % Remove the point from the unconstrained set
VARCnew.POINTSRpw = [VARCnew.POINTSRpw; ipoint_control] ;  % Add it to the weight/position constrained set
% Now we must assign  to xNEW and wNEW the constrained values
% In the case of xNEW, we do not need do to anything (it remains in the
% same position).
% As for the weights, we simply make
wcontr_old = wNEW(ipoint_control) ;
iterWEIGHTS = 1 ;
SALIRloc = 0 ;
while iterWEIGHTS <= DATALOC.SECOND_STAGE_ITERATIONS  && SALIRloc == 0
    % This is the loop in which the weight of the "control" point is
    % gradually diminished     
    wNEW(ipoint_control) = wcontr_old*(1-iterWEIGHTS/DATALOC.SECOND_STAGE_ITERATIONS) ;
    disp('----------------------------------------')
    disp(['ITER = ',num2str(iterWEIGHTS),' of ',num2str(DATALOC.SECOND_STAGE_ITERATIONS),' Reducing weight point ',num2str(ipoint_control),' from ',num2str(wcontr_old),' to ',num2str(wNEW(ipoint_control) ),' (attempt nº ',num2str(iremove),')'])
    disp('---------------------------------------')
         
    [xNEW,wNEW,VARCnew,SALIRloc,SALIR,DATALOC ] =...
        NEWTONRmod(normB,VARCnew,wNEW,DATALOC,TOL,b,xNEW,VAR_SMOOTH_FE,...
        POLYINFO) ;    
    if SALIRloc == 0
        DATALOC.HISTORY_LOCAL.x{end+1} = xNEW ;
        DATALOC.HISTORY_LOCAL.w{end+1}  = wNEW ;
        if iterWEIGHTS == 1
            DATALOC.HISTORY_LOCAL.ControlPoints= ipoint_control ;
        end
        iterWEIGHTS = iterWEIGHTS + 1;       
    end
end
