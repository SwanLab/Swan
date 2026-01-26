function  [rotMATloc,e1] = LocalRotationMatrixNEW(xSIGN,xFIN,xINI,e3,rotX,...
    DATAIN,xMIN,ROT_LOC_SLICE)

% This function should be clean up. It is really a mess !!! JAHO,
% 12-Jan-2019

if  nargin == 0
    load('tmp.mat')
end

   

    
    % Second vector
    if length(xFIN) == 3 && (isempty(DATAIN.ELEVATION_Z) || DATAIN.ELEVATION_Z ==0)  
        e1 = xSIGN*(xFIN-xINI)'/norm(xFIN-xINI) ;
        e2 = cross(e3,e1) ;
        e2 = e2/norm(e2) ;
        % Third vector
        e3 = cross(e1,e2)  ;
        rotMATloc = [e1,e2,e3]*ROT_LOC_SLICE' ;
        % Rotation around x-axis
        rotMATloc= rotMATloc*rotX ;
    else
        % Plane, curved  and helicoidal beams
         e1 = xSIGN*(xFIN(1:2)-xINI(1:2))'/norm(xFIN(1:2)-xINI(1:2)) ;
        e2 = e1(1:2) ;
        e2(1) = -e1(2) ;
        e2(2) = e1(1) ;
        if  length(xFIN) == 2
        rotMATloc = [e1,e2];
       
        else 
            e1 = [e1; 0] ; 
            e2 = [e2; 0] ; 
            e3 = cross(e1,e2) ; 
            rotMATloc = [e1,e2,e3];
        end
        
         rotMATloc =rotMATloc*ROT_LOC_SLICE' ;
    end
    
    % TWISTING
    
    if DATAIN.ISTWIST_ANGLE == 1
        s = (xINI(1) -xMIN)/(xFIN(1)-xINI(1)) ;
        ANGLE_TOTAL = DATAIN.angDOM*s ;
        TWIST_ROTATION = [1 0 0
            0 cos(ANGLE_TOTAL)  -sin(ANGLE_TOTAL)
            0 sin(ANGLE_TOTAL)  cos(ANGLE_TOTAL)] ;
        rotMATloc = rotMATloc*TWIST_ROTATION ;
    end
    