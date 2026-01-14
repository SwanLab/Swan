function  [NAMEPROJECTS,NAME_PROJ_LOC,NUMBER_OF_DOMAINS,IMPOSED_DISP,BoundaryConditions] ...
    = TrainingSet1D(NUMBER_OF_SLICES,NAME_ROOT_FE_BEAM,PRESCRIBED_VALUE,DATARUN)

NAME_PROJ_LOC = {'axial','pbend_z','sbend_z'} ;

NAMEPROJECTS = []; 
NUMBER_OF_DOMAINS = [] ; 
IMPOSED_DISP = [] ; 
BoundaryConditions = [] ; 
if ~isempty(PRESCRIBED_VALUE)
%DATA.INCLUDE_SELFWEIGTH = 1;
%%%%%%%%%%%%%%%%%
nprojects = 3 ;
IMPOSED_DISP.RIGHT_END = zeros(nprojects,3) ;
IMPOSED_DISP.LEFT_END = zeros(nprojects,3) ;

% 
% GENERALIZED_FORCES_ENDS_BEAM.RIGHT_END = zeros(nprojects,6) ;  % Generalized forces
% GENERALIZED_FORCES_ENDS_BEAM.LEFT_END = zeros(nprojects,6) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ndom = NUMBER_OF_SLICES;

BoundaryConditions = cell(size(NAME_PROJ_LOC)) ;
%BoundaryConditions(1:4) = {'PERIODIC_BEAMS_MASS_MATRIX'};

DATARUN = DefaultField(DATARUN,'TYPE_BOUNDARY_CONDITIONS','PERIODIC_BEAMS_MASS_MATRIX') ;

BoundaryConditions(:) = {DATARUN.TYPE_BOUNDARY_CONDITIONS} ;


% %---------------------------------------
% VALUE = 1e3 ; % Value force
% LENGTH_MOMENT = 2 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%
iproj = 1; % Axial test
NAMEPROJECTS{iproj} = [NAME_ROOT_FE_BEAM,NAME_PROJ_LOC{iproj}];
IMPOSED_DISP.RIGHT_END(iproj,1)= PRESCRIBED_VALUE.DISPLACEMENT;
NUMBER_OF_DOMAINS(iproj) = ndom	;

%---------------------------------------
iproj = 2; % Pure bending test -z
NAMEPROJECTS{iproj} =  [NAME_ROOT_FE_BEAM,NAME_PROJ_LOC{iproj}];;
IMPOSED_DISP.RIGHT_END(iproj,3)= PRESCRIBED_VALUE.ROTATION;
IMPOSED_DISP.LEFT_END(iproj,3)= -PRESCRIBED_VALUE.ROTATION;
NUMBER_OF_DOMAINS(iproj) = ndom	;

  
%----%------- Simple bending test z
iproj =3;
NAMEPROJECTS{iproj} =  [NAME_ROOT_FE_BEAM,NAME_PROJ_LOC{iproj}]; ;
IMPOSED_DISP.RIGHT_END(iproj,3)= PRESCRIBED_VALUE.ROTATION;
 
 

NUMBER_OF_DOMAINS(iproj) = ndom	;

else
NAME_PROJ_LOC={} ; % 27-May-2019
end