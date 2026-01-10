function [ROT_LOC_SLICE,rotX,rr1,e3] =  DefinitionZaxis(PROP,MESH1D,CENTRf1,CENTRf2,DATAIN)

% Definition of local z-axis
% ---------------------------
zLOCALdef  = {'PLANE','zGLOBAL',0} ;  % Default value
PROP = DefaultField(PROP,'zLOCAL',zLOCALdef) ; %  Direction local z
% Available options
% MESH1D.PROP(itypeBEAM).zLOCAL{1}  ='PLANE'.  Local z defined by given the unit normal of a given plane. Default option

% MESH1D.PROP(itypeBEAM).zLOCAL{2}  ='zGLOBAL' or itypeBEAM2. If = 'zGLOBAL', the local z
% axis is parallel to the global one.
% If = itypeBEAM2, then it is assumed that elements of type itypeBEAM1 and itypeBEAM2
% lie on a plane. Accordingly, the local z vector will be normal to such a
% plane

% MESH1D.PROP(itypeBEAM).zLOCAL{3} --> Once the z-local axis has been defined by the preceding two arguments, a rotation given
% in degrees by the third argument can be carried out

% DEFAULT VALUES
PROP = DefaultField(PROP,'zLOCAL',{'PLANE','zGLOBAL',0}) ;
if isempty( PROP.zLOCAL )
    PROP.zLOCAL = {'PLANE','zGLOBAL',0} ;
end



TYPE_ROTATION = PROP.zLOCAL{1} ;
ARG2 = PROP.zLOCAL{2} ;
ANGLE = PROP.zLOCAL{3} ;

switch  TYPE_ROTATION
    case 'PLANE'
        if ischar(ARG2)
            if strcmp(ARG2,'zGLOBAL')
                e3 = [0,0,1]' ;
            else
                error('Option not implemented')
            end
        else
            % One element of the type under study
            ielem1 =  INDEXelem(1) ;
            nodoINI = MESH1D.CN(ielem1,1) ;
            nodoFIN = MESH1D.CN(ielem1,2) ;
            xINI = MESH1D.COOR(nodoINI,:) ;
            xFIN = MESH1D.COOR(nodoFIN,:) ;
            vector1 = (xFIN-xINI)/norm(xFIN-xINI) ;
            % One element of the other type
            INDEXelem2 = find(MESH1D.MaterialType == ARG2) ;
            ielem2 = INDEXelem2(2) ;
            nodoINI = MESH1D.CN(ielem2,1) ;
            nodoFIN = MESH1D.CN(ielem2,2) ;
            xINI = MESH1D.COOR(nodoINI,:) ;
            xFIN = MESH1D.COOR(nodoFIN,:) ;
            vector2 = (xFIN-xINI)/norm(xFIN-xINI) ;
            e3 = cross(vector1,vector2) ;
            
            if norm(e3)<=1e-10
                error('Ill-defined plane !! Vectors of both set are aligned')
            end
            e3 = -e3'/norm(e3) ;
            
        end
    otherwise
        error('Option not implemented')
end


%e3 = rotX*e3 ;
% Rotation around x-axis,
% ------------------------
% Default value is  the  identity

rotX = [1    0       0
    0  cosd(ANGLE) -sind(ANGLE)
    0  sind(ANGLE) cosd(ANGLE)      ] ;




% -----------------------------------------------------------------
% Relative rotation  reference configuration
% Helicoidal structures


% if isempty(DATAIN.ELEVATION_Z) || DATAIN.ELEVATION_Z ==0


if length(CENTRf2)==3 && (isempty(DATAIN.ELEVATION_Z) || DATAIN.ELEVATION_Z ==0)
    rr1 = (CENTRf2-CENTRf1)/norm(CENTRf2-CENTRf1) ;
    rr3 = [0 0 1]' ;
    % Second vector  is perpendicular to r1 and nr
    rr2 = cross(rr3,rr1) ;
    % Third vector
    rr3 = cross(rr1,rr2) ;
    ROT_LOC_SLICE = [rr1,rr2,rr3] ; %
else
    rr1 = (CENTRf2(1:2)-CENTRf1(1:2))/norm(CENTRf2(1:2)-CENTRf1(1:2)) ;
    
    rr2 = zeros(2,1) ;
    rr2(1) = -rr1(2);
    rr2(2) = rr1(1)  ;
    if length(CENTRf2)==2
        ROT_LOC_SLICE = [rr1,rr2 ] ; %
    else
        % Helicoidal structure
        rr1 = [rr1; 0] ;
        rr2 = [rr2; 0] ;
        rr3 = cross(rr1,rr2) ;
        ROT_LOC_SLICE = [rr1,rr2,rr3] ; %
    end
end

%