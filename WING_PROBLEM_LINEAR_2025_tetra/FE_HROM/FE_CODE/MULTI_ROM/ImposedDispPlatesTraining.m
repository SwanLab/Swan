function [NAMEPROJECTS,IMPOSED_DISP,NUMBER_OF_DOMAINS] = ....
    ImposedDispPlatesTraining(PRESCRIBED_VALUE,NUMBER_OF_RVES,NAME_ROOT_FE_BEAM,DATARUN)

if nargin == 0
    load('tmp3.mat')
end
% Defining set of FE training tests

% --------------------------------------
% Fixed faces while moving one of them
%-------------------------------------

 

FACES_WITH_NON_ZERO_DISPLACEMENTS = DATARUN.FACES_WITH_NON_ZERO_DISPLACEMENTS ;

 

DATARUN = DefaultField(DATARUN,'PRESCRIBED_DISPLACEMENTS',1) ;

[NAMEPROJECTS,IMPOSED_DISP,NUMBER_OF_DOMAINS] = ....
    LocalNamesPRESCRIBED(PRESCRIBED_VALUE,NUMBER_OF_RVES,NAME_ROOT_FE_BEAM,FACES_WITH_NON_ZERO_DISPLACEMENTS,DATARUN) ;



switch DATARUN.BoundaryConditionType
    case 'MIXED_LINEAR_AND_ZERO_PRESCRIBED'
        
        FACES_WITH_NON_ZERO_DISPLACEMENTS = [3,4] ;
        NAME_ROOT_FE_BEAM = [NAME_ROOT_FE_BEAM,'_2_'] ;
        
        [NAMEPROJECTS_2,IMPOSED_DISP_2,NUMBER_OF_DOMAINS_2] = ....
            LocalNamesPRESCRIBED(PRESCRIBED_VALUE,NUMBER_OF_RVES,NAME_ROOT_FE_BEAM,FACES_WITH_NON_ZERO_DISPLACEMENTS,DATARUN) ;
        
        nproj = length(NAMEPROJECTS)+length(NAMEPROJECTS_2);
        NAMEPROJECTS_new =cell(1,nproj) ;
        NAMEPROJECTS_new(1:length(NAMEPROJECTS)) =  NAMEPROJECTS;
        NAMEPROJECTS_new(length(NAMEPROJECTS)+1:end)  =  NAMEPROJECTS_2 ;
        
        if DATARUN.PRESCRIBED_DISPLACEMENTS == 1
            IMPOSED_DISP_new =cell(size(IMPOSED_DISP)) ;
            for iface = 1:length(IMPOSED_DISP_new)
                IMPOSED_DISP_new{iface} = [IMPOSED_DISP{iface}; IMPOSED_DISP_2{iface}] ;
            end
            IMPOSED_DISP = IMPOSED_DISP_new ;
            NUMBER_OF_DOMAINS_new =cell(1,nproj) ;
            NUMBER_OF_DOMAINS_new(1:length(NAMEPROJECTS)) =  NUMBER_OF_DOMAINS;
            NUMBER_OF_DOMAINS_new(length(NAMEPROJECTS)+1:end)  =  NUMBER_OF_DOMAINS_2 ;
            NUMBER_OF_DOMAINS = NUMBER_OF_DOMAINS_new ;
            NAMEPROJECTS  = NAMEPROJECTS_new ;
        else
            NAMEPROJECTS  = NAMEPROJECTS_new ;
            
        end
        
end

end
