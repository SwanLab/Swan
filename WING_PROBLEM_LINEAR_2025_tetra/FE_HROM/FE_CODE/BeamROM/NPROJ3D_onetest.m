function [NAMEPROJECTS,NAME_PROJ_LOC,NUMBER_OF_DOMAINS,IMPOSED_DISP,BoundaryConditions,GENERALIZED_FORCES_ENDS_BEAM] ...
    = NPROJ3D_onetest(NUMBER_OF_SLICES,NAME_ROOT_FE_BEAM,PRESCRIBED_VALUE,DATARUN)


%DATA.INCLUDE_SELFWEIGTH = 1;
%%%%%%%%%%%%%%%%%

NAME_PROJ_LOC = {'onetest'} ;
NAMEPROJECTS = [] ;
NUMBER_OF_DOMAINS = [] ;
IMPOSED_DISP = [] ;
BoundaryConditions = [] ;
GENERALIZED_FORCES_ENDS_BEAM = [] ; 



if ~isempty(PRESCRIBED_VALUE)
    nprojects = 1 ;
    IMPOSED_DISP.LEFT_END = zeros(nprojects,6) ;
    IMPOSED_DISP.RIGHT_END = zeros(nprojects,6) ;
 %   GENERALIZED_FORCES_ENDS_BEAM{1} = zeros(nprojects,6) ;  % Generalized forces
 %   GENERALIZED_FORCES_ENDS_BEAM{2} = zeros(nprojects,6) ;
    %
    % GENERALIZED_FORCES_ENDS_BEAM.RIGHT_END = zeros(nprojects,6) ;  % Generalized forces
    % GENERALIZED_FORCES_ENDS_BEAM.LEFT_END = zeros(nprojects,6) ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ndom = NUMBER_OF_SLICES;
    
    BoundaryConditions = cell(size(NAME_PROJ_LOC)) ;
    %BoundaryConditions(1:4) = {'PERIODIC_BEAMS_MASS_MATRIX'};
    %DATARUN = DefaultField(DATARUN,'TYPE_BOUNDARY_CONDITIONS','PRESCRIBED_ENDS_BEAMS')  ;
    BoundaryConditions(:) = {'PRESCRIBED_ENDS_BEAMS'} ;
    
    PRESCRIBED_VALUE = DefaultField(PRESCRIBED_VALUE,'GENERALIZED_DISPLACEMENTS',[0,0,0,0,0,0]) ; 

    % %---------------------------------------
    % VALUE = 1e3 ; % Value force
    % LENGTH_MOMENT = 2 ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    iproj = 1; %  
    NAMEPROJECTS{iproj} = [NAME_ROOT_FE_BEAM,NAME_PROJ_LOC{iproj}];
    IMPOSED_DISP.RIGHT_END(iproj,:)= PRESCRIBED_VALUE.GENERALIZED_DISPLACEMENTS;
%     VAL = 1e3 ; 
%     GENERALIZED_FORCES_ENDS_BEAM{2} = PRESCRIBED_VALUE.FORCES   ;   
%    % Generalized forces

    NUMBER_OF_DOMAINS(iproj) = ndom	;
      
    
end