function Dloc = RigidBodyMotionALL(TRANSLATION,INTERVAL,TIMELOC,rnodLOC,Dloc,ndim,ROTATION,COORfaceREL)
%--------------------------------------------------------------------------
% Dloc = RigidBodyMotionALL(TRANSLATION,INTERVAL,TIMELOC,rnodLOC,Dloc,ndim,ROTATION,COORfaceREL)
%
% PURPOSE:
%   Computes the rigid-body displacement vector (translation + rotation)
%   at a given time step for all nodes on a prescribed surface.
%
%   Both translation and rotation are activated only within a given time 
%   interval. The result is accumulated into the current displacement 
%   vector Dloc.
%
% INPUT:
%   TRANSLATION - Struct with translation info (AMPLITUDE, TIMEFUN, etc.).
%   INTERVAL    - [t_start, t_end], activation window for the motion.
%   TIMELOC     - Current time step value.
%   rnodLOC     - List of nodes belonging to the surface.
%   Dloc        - Current displacement vector (to be updated).
%   ndim        - Number of spatial dimensions (2 or 3).
%   ROTATION    - Struct with rotation info (axes, max angles, time laws).
%   COORfaceREL - Coordinates of surface nodes relative to centroid.
%
% OUTPUT:
%   Dloc - Updated displacement vector for the prescribed rigid-body motion.
%
% FUNCTIONALITY:
%   - Adds pure translation displacements to Dloc if TRANSLATION is defined.
%   - Adds rotation-induced displacements:
%       • 2D: calls DispRotation_z_2D (rotation about z-axis).
%       • 3D: calls DispRotation (general 3D rotations).
%
% NOTES:
%   - The function uses helper routines AmplitudeTranslation, DispRotation_z_2D, 
%     and DispRotation to compute motion increments.
%   - Previous “manual” formulas for rotations (commented out) are replaced by 
%     modular helper functions for clarity and efficiency.
%
% JAHO
%--------------------------------------------------------------------------

% --- Translation contribution --------------------------------------------
if ~isempty(TRANSLATION)
    % Compute translation amplitude at this time step (all nodes move equally)
    AMPLITUDE = AmplitudeTranslation(TRANSLATION,INTERVAL,TIMELOC,rnodLOC);
    
    % Add translation displacement to current displacement vector
    Dloc = Dloc + AMPLITUDE;
end

% --- Rotation contribution -----------------------------------------------
if ~isempty(COORfaceREL)
    if ndim == 2
        % 2D case: rotation around the out-of-plane z-axis
        dROTATION = DispRotation_z_2D(ROTATION,INTERVAL,TIMELOC,COORfaceREL);
        Dloc = Dloc + dROTATION(:);
    else
        % 3D case: general rotation (sequential rotations about x,y,z axes)
        dROTATION = DispRotation(ROTATION,INTERVAL,TIMELOC,COORfaceREL);
        Dloc = Dloc + dROTATION(:);
    end
end

% Before being commented by ChatGPT 5, 24th August 2025
% function  Dloc = RigidBodyMotionALL(TRANSLATION,INTERVAL,TIMELOC,rnodLOC,Dloc,ndim,ROTATION,COORfaceREL)
% 
% % TRANSLATION
% 
% if ~isempty(TRANSLATION)
%     AMPLITUDE = AmplitudeTranslation(TRANSLATION,INTERVAL,TIMELOC,rnodLOC) ;
%     %         FactorSteps =  TRANSLATION.TIMEFUN(TIMELOC) ;
%     %         FactorSteps = FactorSteps.*(TIMELOC >= INTERVAL(1)) ;
%     %         FactorSteps = FactorSteps.*(TIMELOC <= INTERVAL(2)) ;
%     %         AMPLITUDE= FactorSteps*TRANSLATION.AMPLITUDE' ;
%     %         AMPLITUDE = repmat(AMPLITUDE,length(rnodLOC),1) ;
%     Dloc = Dloc  + AMPLITUDE ;
% end
% % ROTATION
% if ~isempty(COORfaceREL)
%     if ndim == 2
%         dROTATION = DispRotation_z_2D(ROTATION,INTERVAL,TIMELOC,COORfaceREL) ;
%         %             FactorSteps =  ROTATION.Z.TIMEFUN(TIMELOC) ;
%         %             FactorSteps = FactorSteps.*(TIMELOC >= INTERVAL(1)) ;
%         %             FactorSteps = FactorSteps.*(TIMELOC <= INTERVAL(2)) ;
%         %             ANGLE= FactorSteps*ROTATION.Z.MAXANGLE ;
%         %             R = [cosd(ANGLE) -sind(ANGLE)
%         %                 sind(ANGLE)  cosd(ANGLE)          ] ;
%         %             d =   R*COORfaceREL'  -  COORfaceREL' ;
%         Dloc =Dloc + dROTATION(:) ;
%     else
%         
%         %             R =eye(3) ;
%         %             % Rotations 3D
%         %             ROTAROUND_glo = {'Z','Y','Z'} ;
%         %             IND1_glo = [3,2,1];
%         %             for indrotation =1:3
%         %                 IND1 = IND1_glo(indrotation) ;
%         %                 INDblock = setdiff(1:3,IND1) ;
%         %                 ROTAROUND = ROTAROUND_glo{indrotation} ;
%         %                 FactorSteps =  ROTATION.(ROTAROUND).TIMEFUN(TIMELOC) ;
%         %                 FactorSteps = FactorSteps.*(TIMELOC >= INTERVAL(1)) ;
%         %                 FactorSteps = FactorSteps.*(TIMELOC <= INTERVAL(2)) ;
%         %                 ANGLE= FactorSteps*ROTATION.(ROTAROUND).MAXANGLE ;
%         %                 MAT = [cosd(ANGLE) -sind(ANGLE)
%         %                     sind(ANGLE)  cosd(ANGLE)] ;
%         %                 Rz  =zeros(3) ;
%         %                 Rz(INDblock,INDblock) =   MAT ;
%         %                 Rz(IND1,IND1) =   1 ;
%         %                 R =Rz*R ;
%         %             end
%         %             d =   R*COORfaceREL'  -  COORfaceREL' ;
%         dROTATION = DispRotation(ROTATION,INTERVAL,TIMELOC,COORfaceREL);
%         Dloc =Dloc + dROTATION(:) ;
%         
%     end
%     
% end
% 
%  