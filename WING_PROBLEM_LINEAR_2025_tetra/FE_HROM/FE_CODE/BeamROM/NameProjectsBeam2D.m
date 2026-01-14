function [NAMEPROJECTS,NAME_PROJ_LOC,NUMBER_OF_DOMAINS,IMPOSED_DISP,BoundaryConditions] ...
    = NameProjectsBeam2D(NUMBER_OF_SLICES,NAME_ROOT_FE_BEAM,PRESCRIBED_VALUE,DATARUN)

if nargin == 0
    load('tmp4.mat')
end

DATARUN = DefaultField(DATARUN,'TYPE_TRAINING','NORMAL'); 


switch DATARUN.TYPE_TRAINING 
    case 'NORMAL'
        [NAMEPROJECTS,NAME_PROJ_LOC,NUMBER_OF_DOMAINS,IMPOSED_DISP,BoundaryConditions] ...
    = TrainingSet1D(NUMBER_OF_SLICES,NAME_ROOT_FE_BEAM,PRESCRIBED_VALUE,DATARUN) ; 
    case 'ONLY_AXIAL'
         [NAMEPROJECTS,NAME_PROJ_LOC,NUMBER_OF_DOMAINS,IMPOSED_DISP,BoundaryConditions] ...
    = TrainingSet1DAxial(NUMBER_OF_SLICES,NAME_ROOT_FE_BEAM,PRESCRIBED_VALUE,DATARUN) ; 
case 'ONLY_PURE_BENDING'
         [NAMEPROJECTS,NAME_PROJ_LOC,NUMBER_OF_DOMAINS,IMPOSED_DISP,BoundaryConditions] ...
    = TrainingSet1D_pbend(NUMBER_OF_SLICES,NAME_ROOT_FE_BEAM,PRESCRIBED_VALUE,DATARUN) ;

case 'NONE'
        NAMEPROJECTS = [] ; 
        NAME_PROJ_LOC = [] ; 
        NUMBER_OF_DOMAINS = [] ; 
        IMPOSED_DISP = [] ; 
        BoundaryConditions = [] ; 
    otherwise
        error('Option Not Implemented')
        
    
end
