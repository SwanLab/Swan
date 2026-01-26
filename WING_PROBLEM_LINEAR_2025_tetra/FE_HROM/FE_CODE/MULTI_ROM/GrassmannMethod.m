function  [ERROR_APPROX,INTERPOLATED_BASES ] = GrassmannMethod(SCALE_FACTOR,...
    FACTOR_TEST,NAME_MODES,FACTOR_SCALES_ALL,TYPE,DATAINPUT)

if nargin == 0
    load('tmp.mat')
end
%%%%%%%%%%%55
FOLDER = ['MODES',filesep] ;


% Collecting modes several projects
ModesMatrix = CollectModesProjects(NAME_MODES,FOLDER,TYPE,DATAINPUT) ;

if length(ModesMatrix)==1
    error('The number of matrices shoulb be greater than 2')
elseif length(ModesMatrix)==2
    INTERPOLATED_BASES =  InterpolationTwoSubspaces(ModesMatrix,FACTOR_SCALES_ALL,FACTOR_TEST) ;
    
else
    INTERPOLATED_BASES =  InterpolationSeveralSubspaces...
        (ModesMatrix,FACTOR_SCALES_ALL,FACTOR_TEST,SCALE_FACTOR,DATAINPUT) ;
    
end



 

%%% Comparison with the modes  obtained by running the FE program
% ----------------------------
RUN_MODES.NAME_MODES_FILES = DATAINPUT.NAME_WS_TEST;
FOLDER = [cd,filesep,'MODES',filesep] ;

NAME_MODES   =[FOLDER,'MODES_',RUN_MODES.NAME_MODES_FILES,'.mat'] ;

if ~exist(NAME_MODES ,'file')
    SCALE_FACTOR_loc = FACTOR_TEST  ;
  %  SCALE_FACTOR_loc.X = SCALE_FACTOR.X ;
    
    RUN_FE.SCALE_FACTOR = SCALE_FACTOR_loc ;
    LSCRIPT_fun(DATAINPUT.EXECUTABLE_FOLDER,DATAINPUT.NAME_GEOMETRY_MATERIAL_DATA,...
        DATAINPUT.NAME_LOAD_DATA,RUN_MODES,RUN_FE) ;
end


 switch  TYPE
        case {'BRED'}
            load(NAME_MODES,'BASES','Bdom') ;
            Udisp =  BASES.DISPLACEMENTS.U ; 
            ModesMatrixFE = Bdom*Udisp ;
     otherwise
            load(NAME_MODES,'BASES','DATA_REFMESH') ;

            ModesMatrixFE = BASES.(TYPE).U ;
    end


% ERROR 
ERROR_APPROX = INTERPOLATED_BASES-ModesMatrixFE*(ModesMatrixFE\INTERPOLATED_BASES) ; 
ERROR_APPROX = sum(ERROR_APPROX.^2,1) ; 
normERR = sum(ModesMatrixFE.^2,1) ; 
ERROR_APPROX = sqrt(ERROR_APPROX./normERR)*100 ;

disp(['Approximation error (%) = ',num2str(ERROR_APPROX)])

% 