function [xNEW,wNEW,ISNEGATIVE] = UpdateSolutionWeightsPosition(xNEW,wNEW,delta_q,DATALOC) 


%New point
q_kp1 = [xNEW(:);wNEW] +delta_q ;
m = length(wNEW);
ISNEGATIVE = 0 ;

ndim = size(xNEW,2) ;

wNEW = q_kp1(ndim*m+1:end) ;
if ndim == 2
    xNEW = [q_kp1(1:m) , q_kp1(m+1:2*m) ] ;
    xxLIM =DATALOC.xLIM(1,:) ;
    yyLIM = DATALOC.xLIM(2,:) ;
    
    IS_INSIDE_DOMAIN =   (all(xNEW(:,1) >= xxLIM(1)) && all(xNEW(:,1) <= xxLIM(2)) && ...
        all(xNEW(:,2) >= yyLIM(1)) && all(xNEW(:,2) <= yyLIM(2))) ;
    
    if    IS_INSIDE_DOMAIN
        
        % disp('All points are admissible !!!!!!!')
    else
        %dbstop('34')
        disp('Point out of the domain ...')
        ISNEGATIVE = 1 ;
    end
else
    xNEW = [q_kp1(1:m) , q_kp1(m+1:2*m),  q_kp1(2*m+1:3*m) ] ;
    
    
    xxLIM =DATALOC.xLIM(1,:) ;
    yyLIM = DATALOC.xLIM(2,:) ;
    zzLIM = DATALOC.xLIM(3,:) ;
    
    IS_INSIDE_DOMAIN =  (all(xNEW(:,1) >= xxLIM(1)) && all(xNEW(:,1) <= xxLIM(2)) && ...
        all(xNEW(:,2) >= yyLIM(1)) && all(xNEW(:,2) <= yyLIM(2))  && ...
        all(xNEW(:,3) >= zzLIM(1)) && all(xNEW(:,3) <= zzLIM(2)))  ;
    
    
    if      IS_INSIDE_DOMAIN
        
    else
        %dbstop('34')
        disp('Point out of the domain ...')
        ISNEGATIVE = 1 ;
    end
end

