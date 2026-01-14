function [VARCnew,xNEW,wNEW,wcontr_old,ipoint_control ]= EliminationStrategySTEP(indREM,iremove,VARCnew,xNEW,wNEW,DATALOC,iterWEIGHTS,wcontr_old,ipoint_control)

if nargin == 0
    load('tmp3.mat')
end

%----------------------------------------------------
iremovLOC = indREM(iremove) ; % Indexes of point belonging to  POINTS_F which is constrained now (because we are going to set its weight to zero,

if iterWEIGHTS == 1
    % and the position will remain unchanged)
    POINTS_F = [VARCnew.POINTSl(:); VARCnew.POINTSRp(:)];
    ipoint_control = POINTS_F(iremovLOC) ;  % Index of the constrained point (global)
    % --------------------------------------------------------------------------------------
    % We would have to guess now whether this point belongs to VARCnew.POINTSl or
    % VARCnew.POINTSRp. However, in this strategy, VARCnew.POINTSRp is always
    % empty. Therefore, we simply make
 %   if iremovLOC <= length(VARCnew.POINTSl)
    VARCnew.POINTSl(iremovLOC)=[]  ;  % Remove the point from the unconstrained set
    %POINTS_F=[VARCnew.POINTSl(:)'; VARCnew.POINTSRp(:)' ];  ;  % Remove the point from the unconstrained set
    
    VARCnew.POINTSRpw = [VARCnew.POINTSRpw; ipoint_control] ;  % Add it to the weight/position contrained set
    % Now we must assign  to xNEW and wNEW the constrained values
    % In the case of xNEW, we do not need do to anything (it remains in the
    % same position).
    % As for the weights, we simply make
    wcontr_old = wNEW(ipoint_control) ;
    
end
wNEW(ipoint_control) = wcontr_old*(1-iterWEIGHTS/DATALOC.SECOND_STAGE_ITERATIONS) ;

disp('----------------------------------------')
disp(['ITER = ',num2str(iterWEIGHTS),' of ',num2str(DATALOC.SECOND_STAGE_ITERATIONS),' Reducing weight point ',num2str(ipoint_control),' from ',num2str(wcontr_old),' to ',num2str(wNEW(ipoint_control) ),' (irem =',num2str(iremove),')'])
disp('---------------------------------------')