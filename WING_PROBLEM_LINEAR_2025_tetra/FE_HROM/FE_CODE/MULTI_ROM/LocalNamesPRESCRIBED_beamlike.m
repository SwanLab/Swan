function   [NAMEPROJECTS,IMPOSED_DISP,NUMBER_OF_DOMAINS,NAMELOC] = ....
    LocalNamesPRESCRIBED_beamlike(PRESCRIBED_VALUE,NUMBER_OF_RVES,NAME_ROOT_FE_BEAM,FACES_WITH_NON_ZERO_DISPLACEMENTS,...
    DATARUN)

if nargin == 0
    load('tmp.mat')
end

IMPOSED_DISP = [] ;

NUMBER_OF_DOMAINS = [] ;


%NONZERO_FACES = FACES_WITH_NON_ZERO_DISPLACEMENTS ;

INCLUDE_ALL_DIRECTIONS = 0 ;  

if DATARUN.ndim ==3
    ndirections = 6 ;
    if INCLUDE_ALL_DIRECTIONS == 0 
        % In-plane bending excluded 
       % ndirections = 6; 
        DIRECTIONS = [1,2,3,4]  ; 
    else
        DIRECTIONS = 1:6 ; 
    end
else
    ndirections = 3 ;
    DIRECTIONS = 1:3 ; 
end
% %if DATARUN.PRESCRIBED_DISPLACEMENTS == 1 ;
VALUES_FM = zeros(ndirections,1) ;
if ~isempty(PRESCRIBED_VALUE)
    if ndirections ==6
        VALUES_FM(1:3) = PRESCRIBED_VALUE.DISPLACEMENT ;
        VALUES_FM(4:6) = PRESCRIBED_VALUE.ROTATION ;
    else
        VALUES_FM(1:2) = PRESCRIBED_VALUE.DISPLACEMENT ;
        VALUES_FM(3) = PRESCRIBED_VALUE.ROTATION ;
    end
else
    VALUES_FM = [] ;
end
% 
% % EMPTY VALUES 
% VALUES_FM = zeros(ndirections,1) ;


%end
nprojects = length(DIRECTIONS)*2;  % 6 along x, 6 along y
NAMEPROJECTS = cell(1,nprojects) ;
NAMELOC = NAMEPROJECTS ;
%if DATARUN.PRESCRIBED_DISPLACEMENTS == 1 ;
NUMBER_OF_DOMAINS = cell(1,nprojects) ;
NUMBER_OF_DOMAINS(:) = {NUMBER_OF_RVES} ;
%end
iacum = 1 ;
nfaces = 4 ;
%if  DATARUN.PRESCRIBED_DISPLACEMENTS == 1
IMPOSED_DISP = cell(nfaces,1) ;
IMPOSED_DISP(:) = {cell(nprojects,ndirections)} ;

% PRESCRIBED DISPLACEMENTS FACES 1 AND 3 
% --------------------------------------
ifaceFIXglo = [1,2] ; 
ifaceGLO = [3,4] ; 



for iloc = 1:length(ifaceFIXglo)
    ifaceFIX = ifaceFIXglo(iloc) ;
    iface = ifaceGLO(iloc) ;
    
    for idirectionLOC = 1:length(DIRECTIONS)
        idirection = DIRECTIONS(idirectionLOC) ; 
        NAMELOC{iacum} = ['Z',num2str(ifaceFIX),'_f',num2str(iface),'_','dir_',num2str(idirection)] ;
        NAMEPROJECTS{iacum} = [NAME_ROOT_FE_BEAM,'_',NAMELOC{iacum}] ;
        IMPOSED_DISP{iface}(iacum,:) = {0};
        IMPOSED_DISP{ifaceFIX}(iacum,:) ={0} ;
        if ~isempty(VALUES_FM)
         IMPOSED_DISP{iface}{iacum,idirection} = VALUES_FM(idirection) ;
        end
         iacum = iacum  + 1;
    end
end
 

%end

% 

% 
% 


end
