function  [BASES,BdomRED] = GrassmannMethodOperators(SCALE_FACTOR,...
    FACTOR_TEST,NAME_MODES,FACTOR_SCALES_ALL,TYPE,DATAINPUT)

if nargin == 0
    load('tmp.mat')
end
%%%%%%%%%%%55
FOLDER = ['MODES',filesep] ;


%  Basis matrices 
% ----------------
TYPE ={'DISPLACEMENTS','REACTIONS','STRESSES'} ;
BASES = [] ; 
for i=1:length(TYPE)    
    [ERROR_APPROXIMATION,BASES_LOC ]= GrassmannMethod(SCALE_FACTOR,...
        FACTOR_TEST,NAME_MODES,FACTOR_SCALES_ALL,TYPE{i},DATAINPUT) ;    
    BASES.(TYPE{i}) = [];
    BASES.(TYPE{i}).U =BASES_LOC;
    BASES.(TYPE{i}).S =ones(size(BASES_LOC,2),1);  
end

% REduced Strain-Displacement matrix
% ------------------------------------
TYPE ='BRED' ;
[ERROR_APPROXIMATION,BdomRED ]= GrassmannMethod(SCALE_FACTOR,...
    FACTOR_TEST,NAME_MODES,FACTOR_SCALES_ALL,TYPE,DATAINPUT) ;

% Integration Weights

Winterp = InterpolationMatricesDirect(FACTOR_SCALES_ALL,W,FACTOR_TEST,SCALE_FACTOR,DATAINPUT) ; 


 