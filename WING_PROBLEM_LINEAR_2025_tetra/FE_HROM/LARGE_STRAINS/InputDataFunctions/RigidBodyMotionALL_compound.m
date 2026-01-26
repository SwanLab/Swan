function  [Dglo,R] = RigidBodyMotionALL_compound(TRANSLATION,INTERVAL,TIMELOC,rnodLOC,Dloc,ndim,ROTATION,COORfaceREL,dLOCAL)

% TRANSLATION
Dglo = zeros(size(Dloc)) ; 
if ~isempty(TRANSLATION)
    AMPLITUDE = AmplitudeTranslation(TRANSLATION,INTERVAL,TIMELOC,rnodLOC) ;
    %         FactorSteps =  TRANSLATION.TIMEFUN(TIMELOC) ;
    %         FactorSteps = FactorSteps.*(TIMELOC >= INTERVAL(1)) ;
    %         FactorSteps = FactorSteps.*(TIMELOC <= INTERVAL(2)) ;
    %         AMPLITUDE= FactorSteps*TRANSLATION.AMPLITUDE' ;
    %         AMPLITUDE = repmat(AMPLITUDE,length(rnodLOC),1) ;
    Dglo =  AMPLITUDE ;
end
% ROTATION
if ~isempty(COORfaceREL)
    if ndim == 2
        [dROTATION,R] = DispRotation_z_2D(ROTATION,INTERVAL,TIMELOC,COORfaceREL) ;
        %             FactorSteps =  ROTATION.Z.TIMEFUN(TIMELOC) ;
        %             FactorSteps = FactorSteps.*(TIMELOC >= INTERVAL(1)) ;
        %             FactorSteps = FactorSteps.*(TIMELOC <= INTERVAL(2)) ;
        %             ANGLE= FactorSteps*ROTATION.Z.MAXANGLE ;
        %             R = [cosd(ANGLE) -sind(ANGLE)
        %                 sind(ANGLE)  cosd(ANGLE)          ] ;
        %             d =   R*COORfaceREL'  -  COORfaceREL' ;
        Dglo =Dglo + dROTATION(:) ;
        
        d_DEF = R*dLOCAL;
        Dglo =Dglo + d_DEF(:) ;
    else
        
        %             R =eye(3) ;
        %             % Rotations 3D
        %             ROTAROUND_glo = {'Z','Y','Z'} ;
        %             IND1_glo = [3,2,1];
        %             for indrotation =1:3
        %                 IND1 = IND1_glo(indrotation) ;
        %                 INDblock = setdiff(1:3,IND1) ;
        %                 ROTAROUND = ROTAROUND_glo{indrotation} ;
        %                 FactorSteps =  ROTATION.(ROTAROUND).TIMEFUN(TIMELOC) ;
        %                 FactorSteps = FactorSteps.*(TIMELOC >= INTERVAL(1)) ;
        %                 FactorSteps = FactorSteps.*(TIMELOC <= INTERVAL(2)) ;
        %                 ANGLE= FactorSteps*ROTATION.(ROTAROUND).MAXANGLE ;
        %                 MAT = [cosd(ANGLE) -sind(ANGLE)
        %                     sind(ANGLE)  cosd(ANGLE)] ;
        %                 Rz  =zeros(3) ;
        %                 Rz(INDblock,INDblock) =   MAT ;
        %                 Rz(IND1,IND1) =   1 ;
        %                 R =Rz*R ;
        %             end
        %             d =   R*COORfaceREL'  -  COORfaceREL' ;
        [dROTATION,R ]= DispRotation(ROTATION,INTERVAL,TIMELOC,COORfaceREL);
        Dglo =Dglo + dROTATION(:) ;
        
        d_DEF = R*dLOCAL;
        Dglo =Dglo + d_DEF(:) ;
        
    end
    
else
    Dglo =Dglo +  dLOCAL(:); 
    
end

 