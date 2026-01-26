function [VARCnew,xNEW,wNEW ]= EliminationStrategyFIXED(indREM,iremove,VARCnew,xNEW,wNEW)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/07_BEAMS_NOGRAD/ControlledVariables3D.mlx
if nargin == 0
    load('tmp1.mat')
end

%----------------------------------------------------

iremovLOC = indREM(iremove) ; % Indexes of point belonging to  POINTS_F which is constrained now (because we are going to set its weight to zero,
% and the position will remain unchanged)
POINTS_F = [VARCnew.POINTSl(:); VARCnew.POINTSRp(:)];
ipoint_control = POINTS_F(iremovLOC) ;  % Index of the constrained point (global)
% --------------------------------------------------------------------------------------
% We would have to guess now whether this point belongs to VARCnew.POINTSl or
% VARCnew.POINTSRp. 
if iremovLOC <= length(VARCnew.POINTSl)
    % We remove the index from this set 
    VARCnew.POINTSl(iremovLOC)=[]  ;  % Remove the point from the unconstrained set
else
    % It belongs to VARCnew.POINTSp; removing it ...
    iremovLOC_p = iremovLOC-length(VARCnew.POINTSl) ; 
    VARCnew.POINTSRp(iremovLOC_p)=[] ;
end

VARCnew.POINTSRpw = [VARCnew.POINTSRpw; ipoint_control] ;  % Add it to the weight/position contrained set 
% Now we must assign  to xNEW and wNEW the constrained values 
% In the case of xNEW, we do not need do to anything (it remains in the
% same position). 
% As for the weights, we simply make
wNEW(ipoint_control) = 0 ;  

disp('----------------------------------------')
disp(['Removing point =',num2str(ipoint_control),' (irem =',num2str(iremove),')'])
disp('---------------------------------------')