
function [dROTATION,R ]= DispRotation_z_2D(ROTATION,INTERVAL,TIMELOC,COORfaceREL)

FactorSteps =  ROTATION.Z.TIMEFUN(TIMELOC) ;
FactorSteps = FactorSteps.*(TIMELOC >= INTERVAL(1)) ;
FactorSteps = FactorSteps.*(TIMELOC <= INTERVAL(2)) ;
ANGLE= FactorSteps*ROTATION.Z.MAXANGLE ;
R = [cosd(ANGLE) -sind(ANGLE)
    sind(ANGLE)  cosd(ANGLE)          ] ;
dROTATION =   R*COORfaceREL'  -  COORfaceREL' ;


end

