function   [NAMEPROJECTS,IMPOSED_DISP,NUMBER_OF_DOMAINS,NAMELOC,IMPOSED_FORCE] = ....
    LocalNamesPRESCRIBED_forces(PRESCRIBED_VALUE,NUMBER_OF_RVES,NAME_ROOT_FE_BEAM,FACES_WITH_NON_ZERO_DISPLACEMENTS,...
    DATARUN)

if nargin == 0
    load('tmp3.mat')
end


IMPOSED_DISP = [] ;

NUMBER_OF_DOMAINS = [] ;


NONZERO_FACES = FACES_WITH_NON_ZERO_DISPLACEMENTS ;
if DATARUN.ndim ==3
    ndirections = 6 ;
else
    ndirections = 3 ;
end
if DATARUN.PRESCRIBED_DISPLACEMENTS == 1 ;
    VALUES_FM = zeros(ndirections,1) ;
    if ndirections ==6
        VALUES_FM(1:3) = PRESCRIBED_VALUE.DISPLACEMENT ;
        VALUES_FM(4:6) = PRESCRIBED_VALUE.ROTATION ;
    else
        VALUES_FM(1:2) = PRESCRIBED_VALUE.DISPLACEMENT ;
        VALUES_FM(3) = PRESCRIBED_VALUE.ROTATION ;
    end
    
end
nprojects = ndirections*length(NONZERO_FACES) + 1;
NAMEPROJECTS = cell(1,nprojects) ;
NAMELOC = NAMEPROJECTS ;
if DATARUN.PRESCRIBED_DISPLACEMENTS == 1 ;
    NUMBER_OF_DOMAINS = cell(1,nprojects) ;
    NUMBER_OF_DOMAINS(:) = {NUMBER_OF_RVES} ;
end
iacum = 1 ;
nfaces = 4 ;
if  DATARUN.PRESCRIBED_DISPLACEMENTS == 1
    IMPOSED_DISP = cell(nfaces,1) ;
    IMPOSED_DISP(:) = {zeros(nprojects,ndirections)} ;
end
for ifaceLOC = 1:length(NONZERO_FACES)
    iface = NONZERO_FACES(ifaceLOC) ;
    for idirection = 1:ndirections
        NAMELOC{iacum} = ['f',num2str(iface),'_','dir_',num2str(idirection)] ;
        NAMEPROJECTS{iacum} = [NAME_ROOT_FE_BEAM,'_',NAMELOC{iacum}] ;
        if DATARUN.PRESCRIBED_DISPLACEMENTS == 1
            IMPOSED_DISP{iface}(iacum,idirection) = VALUES_FM(idirection) ;
        end
        iacum = iacum  + 1;
    end
end

% Applied force
NAMELOC{iacum} = ['_normal_force'] ;
NAMEPROJECTS{iacum} = [NAME_ROOT_FE_BEAM,'_',NAMELOC{iacum}] ;
try
IMPOSED_FORCE = cell(size(NAMEPROJECTS)) ; 
IMPOSED_FORCE{iacum} = PRESCRIBED_VALUE.FORCE ; 
catch
    IMPOSED_FORCE = [] ; 
end

end
