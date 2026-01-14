function ROTMATRIX = LocalRotMatrix(iface,DATA,DOMS,ndim)

DATA = DefaultField(DATA,'ANGrelGLO',[]) ;
if  isempty(DATA.ANGrelGLO)
    %-----------------------------------------------------------------------------
    if iface==2 && ~isempty(DATA.angROTATION_FACE)
        if length(DATA.angROTATION_FACE ) == 1
            total_angle = DOMS(iface)*DATA.angROTATION_FACE ;
        else
            % Variable curvature
            total_angle = sum(DATA.angROTATION_FACE) ;
        end
        DATA = DefaultField(DATA,'ISTWIST_ANGLE',0) ; % == 0
        if DATA.ISTWIST_ANGLE == 0
            ROTMATRIX = [cos(total_angle),sin(total_angle)
                -sin(total_angle)  ,   cos(total_angle)] ;
            if ndim == 3
                ROTMATRIX_loc = eye(3) ;
                ROTMATRIX_loc(1:2,1:2) = ROTMATRIX ;
                ROTMATRIX = ROTMATRIX_loc ;
            end
        else
            ROTMATRIX = [1,0,0 ;
                0, cos(total_angle),-sin(total_angle)
                0 , sin(total_angle)  ,   cos(total_angle)] ;
        end
        %         uBARloc = ROTMATRIX*reshape(uBARloc,ndim,[]) ;
        %         uBARloc = uBARloc(:) ;
    else
        ROTMATRIX = [];
    end
    %---------------------------------------------------------------------------
else
    % Data generated with % CURVED_BEAMS/ParametricCurvature3D.m
    if iface == 1
        ROTMATRIX  = DATA.ANGiniGLO{1} ;
    elseif iface ==2
        ROTMATRIX  = DATA.ANGiniGLO{end}*DATA.ANGrelGLO{end} ;
    else
        error('Option not valid')
    end
    
    
end