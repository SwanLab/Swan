function   [NAMEPROJECTS,IMPOSED_DISP,NUMBER_OF_DOMAINS,NAMES] = ....
    LocalNamesPRESCRIBED_plates(PRESCRIBED_VALUE,NUMBER_OF_RVES,NAME_ROOT_FE_BEAM,FACES_WITH_NON_ZERO_DISPLACEMENTS,...
    DATARUN)

if nargin == 0
    load('tmp1.mat')
end
% See PRESCRIBED_DISP_PLATES_RVE_fun.m . 
NAMES{1} = 'StrainXX' ;  %Trans x, face 3
NAMES{2} = 'StrainXY' ;  % Trans y, face 3
NAMES{3} = 'StrainXZ' ;  % Trans z, face 3
NAMES{4} = 'TorsionX' ;  % Rot x, face 3
NAMES{5} = 'BendingXZ' ;  % Rot y, face 3
NAMES{6} = 'BendingXY' ;  % Rot z, face 3
%%%
NAMES{7} = 'StrainY' ;  % Trans y, face 4
NAMES{8} = 'StrainYZ' ; % Trans z, face 4
NAMES{9} = 'BendingYZ' ; % Rot X, face 4
NAMES{10} = 'TorsionY' ; % Rot Y, face 4
NAMES{11} = 'BendingYX' ; % Roz Z, face 4
%%%%%
NAMES{12} = 'Null1' ;  % Rot x, face 1
%
NAMES{13} = 'Null2';
NAMES{14}  ='Null3' ;


for iproj = 1:length(NAMES)
        NAMEPROJECTS{iproj} = [NAME_ROOT_FE_BEAM,'_',NAMES{iproj}] ;
      
end







IMPOSED_DISP = [] ;

 



% NONZERO_FACES = FACES_WITH_NON_ZERO_DISPLACEMENTS ;
% ndirections = 6 ;
% if DATARUN.PRESCRIBED_DISPLACEMENTS == 1 ;
%     VALUES_FM = zeros(ndirections,1) ;
%     VALUES_FM(1:3) = PRESCRIBED_VALUE.DISPLACEMENT ;
%     VALUES_FM(4:6) = PRESCRIBED_VALUE.ROTATION ;
%
% end
nprojects = length(NAMES);
NUMBER_OF_DOMAINS = cell(1,nprojects) ;
NUMBER_OF_DOMAINS(:) = {NUMBER_OF_RVES} ;

 



end