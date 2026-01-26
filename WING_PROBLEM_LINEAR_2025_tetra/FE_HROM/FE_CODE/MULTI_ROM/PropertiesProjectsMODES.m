function [DATA_TRAINING,DATAIN ]= PropertiesProjectsMODES(DATAIN,NAMEPROJECTS,FOLDER,COMPUTE_MODES_AGAIN,TOL_SINGULAR_VALUES_Hqr)

RVEs_SEL = DATAIN.RVES_SEL ;
IS_RVE_MODE = ones(size(NAMEPROJECTS)) ;
DATAIN = DefaultField(DATAIN,'nU',cell(size(NAMEPROJECTS)));
DATAIN = DefaultField(DATAIN,'nR',cell(size(NAMEPROJECTS))) ;
for iproj = 1:length(NAMEPROJECTS)
    nU{iproj}= DATAIN.nU{iproj};
    nR{iproj}= DATAIN.nR{iproj};;
    RVEs{iproj} =RVEs_SEL ;
end
TOTAL_NUMBER_OF_MODES = [] ;
% If SLICE{i} is empty, then all slices are included in the SVD
DATA_TRAINING = cell(size(NAMEPROJECTS)) ;  % % Full Path Names .mat files containing FE data (Stiffness matrix...)
for itraj = 1:length(DATA_TRAINING)
    %  These paths are specified when running the FE analyses (FE_TRAIN_fast.m, in this case)
    DATA_TRAINING{itraj} = [ FOLDER,filesep,'DATAFE_TRAIN',filesep,NAMEPROJECTS{itraj},'.mat'] ; %
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of modes of each project for deformational displacements,
DATAIN.NMODES_PROJECT_DISPLACEMENTS =  nU;
% Number of reaction modes
DATAIN.NMODES_PROJECT_REACTIONS =   nR;
DATAIN.DOMAINS_TO_INCLUDE_TRAINING = RVEs;
DATAIN.NMODES_TRUNCATE = TOTAL_NUMBER_OF_MODES ;
DATAIN.IS_RVE_MODE = IS_RVE_MODE ;
DATAIN.COMPUTE_MODES_AGAIN = COMPUTE_MODES_AGAIN ;
DATAIN.TOL_SINGULAR_VALUES_Hqr = TOL_SINGULAR_VALUES_Hqr ;





%%% ADDITIONAL NON-RVE PROJECTS
%%% -----------------------------
% DISABLED, 11-JUNIO-2019
% if ~isempty(NonRVEProjects)
%     error('Program this part !!!! ')
%     
%     NAMEWSloc =    NonRVEProjects.NAMEWS  ;   % Absolute path
%     nPROJADD = length(NAMEWSloc)
%     NonRVEProjects  =DefaultField(NonRVEProjects,'nU',cell(nPROJADD,1)) ;
%     NonRVEProjects  =DefaultField(NonRVEProjects,'nR',cell(nPROJADD,1)) ;
%     NonRVEProjects  =DefaultField(NonRVEProjects,'RVEs',cell(nPROJADD,1)) ;
%     
%     for iaddproj = 1:length(NAMEWSloc)
%         DATA_TRAINING{end+1} = NAMEWSloc{iaddproj} ;
%         nU{end+1} = NonRVEProjects.nU{iaddproj}  ;
%         nR{end+1} = NonRVEProjects.nR{iaddproj}  ;  ;
%         RVEs{end+1} =NonRVEProjects.RVEs{iaddproj}  ;  ;   ;
%         IS_RVE_MODE(end+1) = 0 ;
%     end
%     
%     
%     
%     
% end
