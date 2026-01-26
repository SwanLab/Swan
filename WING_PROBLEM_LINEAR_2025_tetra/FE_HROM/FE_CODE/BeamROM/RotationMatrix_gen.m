function [rotMATloc e3]= RotationMatrix_gen(xINI,xFIN,ANG,e3)

DIST = norm(xFIN-xINI) ;
% Rotation matrix
r1 = (xFIN-xINI)'/DIST ;  % Local x vector expressed in the global axis
% Vector perpendicular to the plane defined by r1 and the global x axis
%

if  ~isempty(e3)
    r2 = cross(e3,r1) ; %
    r2 = r2/norm(r2);
    % Third vector
    r3 = cross(r1,r2) ;
    rotMATloc = [r1,r2,r3] ;
else
    if norm(r1-[1 0 0]') >0
        nr = cross(r1,[1 0 0]') ;
        % Second vector  is perpendicular to r1 and nr
        r2 = cross(nr,r1) ; %
        r2 = r2/norm(r2);
        % Third vector
        r3 = cross(r1,r2) ;
        rotMATloc = [r1,r2,r3] ;
        e3 = [r3] ;
        
        %% Rotation around local x-axis
        if ANG ~=0
            rotX = [1    0       0
                0  cosd(ANG) -sind(ANG)
                0  sind(ANG) cosd(ANG)      ] ;
            rotMATloc = rotX*rotMATloc ;
        end
        
        
    else
        rotMATloc = eye(3) ;
        e3 = [0,0,1]' ;
    end
end
