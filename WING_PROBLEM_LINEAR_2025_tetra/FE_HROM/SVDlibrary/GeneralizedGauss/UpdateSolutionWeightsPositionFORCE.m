function [xNEW,wNEW,ISNEGATIVE] = UpdateSolutionWeightsPositionFORCE(xNEW,wNEW,delta_q,DATALOC)
% See comments in UpdateSolutionWeightsPositionFORCE_aux.mlx
if nargin ==0
    load('tmp2.mat')
end

%New point
xOLD = xNEW ;
q_kp1 = [xNEW(:);wNEW] +delta_q ;
m = length(wNEW);
ISNEGATIVE = 0 ;

ndim = size(xNEW,2) ;
wNEW = q_kp1(ndim*m+1:end) ;

if ndim == 2
    xNEW = [q_kp1(1:m) , q_kp1(m+1:2*m) ] ;
elseif ndim == 3
    xNEW = [q_kp1(1:m) , q_kp1(m+1:2*m),  q_kp1(2*m+1:3*m) ] ;
    %  zzLIM = DATALOC.xLIM(3,:) ;
elseif ndim == 1 
    xNEW = [q_kp1(1:m)  ] ;
else 
    
    error('Option not implemented')
    
end



POINTS_INSIDE_DOMAIN = ones(m,1) ;

for idim = 1:ndim
    POINTS_INSIDE_DOMAIN = POINTS_INSIDE_DOMAIN.*(xNEW(:,idim) >= DATALOC.xLIM(idim,1)).*(xNEW(:,idim) <= DATALOC.xLIM(idim,2)) ;
end

xLIM_EXPANDED = zeros(size(DATALOC.xLIM)) ;

GAP_OVER_1 = DATALOC.TOLERANCE_FORCE_POINTS_TO_REMAIN_WITHIN_THE_DOMAIN;

xLIM_EXPANDED(:,1)= DATALOC.xLIM(:,1) -  GAP_OVER_1*(DATALOC.xLIM(:,2)-DATALOC.xLIM(:,1)) ;
xLIM_EXPANDED(:,2)= DATALOC.xLIM(:,2) +  GAP_OVER_1*(DATALOC.xLIM(:,2)-DATALOC.xLIM(:,1)) ;

ADMISSIBLE_POINTS = ones(m,1) ;

for idim = 1:ndim
    ADMISSIBLE_POINTS = ADMISSIBLE_POINTS.*(xNEW(:,idim) >= xLIM_EXPANDED(idim,1)).*(xNEW(:,idim) <= xLIM_EXPANDED(idim,2)) ;
end

%  if any(wNEW <-0.1*sum(wNEW))
%             disp('Negative weights ')
%  end




if     all(ADMISSIBLE_POINTS)
    % If all points are "admissible points", we check which points are
    % actually within the domain; those which are not are return to their
    % original positions
    IND_POINTS_OUTSIDE = find(POINTS_INSIDE_DOMAIN ==0) ;
    
    if ~isempty(IND_POINTS_OUTSIDE)
        disp('Returning points outside the domain ..')
    xNEW(IND_POINTS_OUTSIDE,:) = xOLD(IND_POINTS_OUTSIDE,:) ; 
    end
    
%     for ipoints = 1:length(IND_POINTS_OUTSIDE)
%         xOLDloc  = xOLD(IND_POINTS_OUTSIDE(ipoints),:) ;
%         xNEWloc =  xNEW(IND_POINTS_OUTSIDE(ipoints),:) ;
%         disp(['Moving point: ',num2str(xNEWloc)])
%         disp(['to point    :',num2str(xOLDloc)])
%         disp('--------------')
%         xNEW(IND_POINTS_OUTSIDE(ipoints),:) =  xOLDloc ;
%     end
    
else
    %dbstop('34')
    disp('Point out of the domain ...')
    ISNEGATIVE = 1 ;
end


