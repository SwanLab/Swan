function [NAMEPROJECTS] = NameProjects_only(NAME_ROOT_FE_BEAM,DATARUN)

% Defining set of FE training tests 
NONZERO_FACES =DATARUN.FACES_WITH_NON_ZERO_DISPLACEMENTS ; 
ndirections = 6 ; 
nprojects = ndirections*length(NONZERO_FACES); 
NAMEPROJECTS = cell(1,nprojects) ;
iacum = 1 ; 
nfaces = 4 ; 
for ifaceLOC = 1:length(NONZERO_FACES) 
    iface = NONZERO_FACES(ifaceLOC) ; 
    for idirection = 1:ndirections
        NAMEPROJECTS{iacum} = [NAME_ROOT_FE_BEAM,'_f',num2str(iface),'_','dir_',num2str(idirection)] ; 
        iacum = iacum  + 1; 
    end
end