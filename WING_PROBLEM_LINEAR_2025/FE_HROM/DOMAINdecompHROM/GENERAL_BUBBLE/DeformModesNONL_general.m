function [BasisUdeform,Mdom,PhiRB,MESH,DATA,PsiRBf,Mintf,DATAoffline] ...
    = DeformModesNONL_general(OPERFE,MESH,DATA,...
    SNAPdisp,DATAcommon,DATAoffline,SNAPdisp_AUX)
% DETERMINATION OF deformational  MODES FOR A GIVEN SUBDOMAIN
% General method
%% JAHO, 3-APR-2024, Sandwichez, C/Mallorca Barcelona

if nargin == 0
    load('tmp1.mat')
end


DATAcommon = DefaultField(DATAcommon,'METHOD_TRAINING_ELASTIC_RANGE',[]) ;


if  ~isempty(DATAcommon.METHOD_TRAINING_ELASTIC_RANGE)
    % ------------------------------------------------
    % Specific method for determining invariant modes + complementary
    % modes, 
    % Apr-3-2024
    % see /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/02_BEAMQ8_IprofE.mlx
    % -----------------------------------------------
    switch  DATAcommon.METHOD_TRAINING_ELASTIC_RANGE
        case  'INVARIANT_MODES_PLUS_COMPLEMENTARY_MODES'
              [BasisUdeform,Mdom,PhiRB,MESH,DATA,PsiRBf,Mintf,DATAoffline]  =...
            DeformModesNONL_invCOMPALL(OPERFE,MESH,DATA,...
            SNAPdisp,DATAcommon,DATAoffline,SNAPdisp_AUX) ;
            
        otherwise
            error('Method not implemented')
        
    end
    
else
        % Methods implemented BEFORE 3th April 2024

    
    if isempty(DATAoffline.METHOD_INVARIANT_DEF_MODES)
        [BasisUdeform,Mdom,PhiRB,MESH,DATA,PsiRBf,Mintf,DATAoffline]  =...
            DeformModesNONL(OPERFE,MESHdom,DATA,...
            SNAPdisp,DATAcommon,DATAoffline,SNAPdisp_AUX) ;
        
    else
        DATAoffline.METHOD_INVARIANT_DEF_MODES = DefaultField(DATAoffline.METHOD_INVARIANT_DEF_MODES,'STAGE_TRAINING',1);
        if DATAoffline.METHOD_INVARIANT_DEF_MODES.STAGE_TRAINING ==1
            % This is the first stage of training when invariant modes are to
            % be determined
            %
            
            DeformModesNONL_invMODES(OPERFE,MESH,DATA,...
                SNAPdisp,DATAcommon,DATAoffline,SNAPdisp_AUX) ;
            
            
            BasisUdeform = [] ;,
            Mdom = [] ;
            PhiRB = [] ;
            PsiRBf =[] ; ,Mintf=[] ;
            
            return
            
        else
            [BasisUdeform,Mdom,PhiRB,MESH,DATA,PsiRBf,Mintf,DATAoffline]  =...
                DeformModesNONL(OPERFE,MESH,DATA,...
                SNAPdisp,DATAcommon,DATAoffline,SNAPdisp_AUX) ;
        end
    end
    
end