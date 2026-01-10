function [dROTATION,R] = DispRotation(ROTATION,INTERVAL,TIMELOC,COORfaceREL)

if nargin == 0
    load('tmp1.mat')
end


R =eye(3) ;

   ROTATION = DefaultField(ROTATION,'RotationVectorGlobal',[]) ; 

if ~isempty(ROTATION.RotationVectorGlobal)
    % The data has been provided in the form of a global 
 
        
        
        FactorSteps =  ROTATION.TIMEFUN(TIMELOC) ;
        FactorSteps = FactorSteps.*(TIMELOC >= INTERVAL(1)) ;
        FactorSteps = FactorSteps.*(TIMELOC <= INTERVAL(2)) ;
         
        MAXANGLE = norm(ROTATION.RotationVectorGlobal) ;  
        % In degrees   
        ANGLE= FactorSteps*MAXANGLE;
        
        % Rodrigues formula 
        if any(ROTATION.RotationVectorGlobal)
         v = ROTATION.RotationVectorGlobal/MAXANGLE ; 
        SPIN = [0 -v(3) v(2)
                v(3) 0  -v(1)
                -v(2) v(1) 0] ; 
            
          R = eye(3) + sind(ANGLE)*SPIN + (1-cosd(ANGLE))*(SPIN*SPIN) ;  
            
        end
    
else
    % Rotations 3D
    ROTAROUND_glo = {'Z','Y','X'} ;
    IND1_glo = [3,2,1];
    for indrotation =1:3
        IND1 = IND1_glo(indrotation) ;
        INDblock = setdiff(1:3,IND1) ;
        ROTAROUND = ROTAROUND_glo{indrotation} ;
        FactorSteps =  ROTATION.(ROTAROUND).TIMEFUN(TIMELOC) ;
        FactorSteps = FactorSteps.*(TIMELOC >= INTERVAL(1)) ;
        FactorSteps = FactorSteps.*(TIMELOC <= INTERVAL(2)) ;
        ANGLE= FactorSteps*ROTATION.(ROTAROUND).MAXANGLE ;
        MAT = [cosd(ANGLE) -sind(ANGLE)
            sind(ANGLE)  cosd(ANGLE)] ;
        Rz  =zeros(3) ;
        Rz(INDblock,INDblock) =   MAT ;
        Rz(IND1,IND1) =   1 ;
        R =Rz*R ;
    end
end

dROTATION =   R*COORfaceREL'  -  COORfaceREL' ;
