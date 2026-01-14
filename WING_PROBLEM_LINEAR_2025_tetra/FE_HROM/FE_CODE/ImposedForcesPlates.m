function FUNinput = ImposedForcesPlates(FUNinput,IMPOSED_FORCE,ndom,NUMBER_OF_DOMAINS) 
if nargin == 0
    load('tmp3.mat')
    IMPOSED_FORCE.RVES_APPLIED  ={[1,2],[1,2]}; 
end

FUNinput.INPUTS.PLATELOADS.LOAD_UNIFORM(:) = {zeros(ndom,3) };
FUNinput.INPUTS.PLATELOADS.ISLOCAL(:) = {zeros(ndom,1) };
%             PRESCRIBED_VALUE.FORCE.VALUE = 0.1 ;   % Prescribed force normal to domain (per unit surface)
% PRESCRIBED_VALUE.FORCE.SURFACE = 5 ;   % Surface on which the force act
% PRESCRIBED_VALUE.FORCE.RVES_APPLIED = {2,2} ;   % RVEs on which force is applied

iSURFACE = IMPOSED_FORCE.SURFACE         ;
VALUE_FORCE = IMPOSED_FORCE.VALUE        ;
RVE_applied = IMPOSED_FORCE.RVES_APPLIED ;
ndomLOCAL = NUMBER_OF_DOMAINS            ;
all_locals = 1:prod(ndomLOCAL) ;
DOMAINS = reshape(all_locals,ndomLOCAL(1),[]); 
DOMAINS_SELECT  = DOMAINS(RVE_applied{1},RVE_applied{2}) ;  
 


FUNinput.INPUTS.PLATELOADS.LOAD_UNIFORM{iSURFACE}(DOMAINS_SELECT,1) = VALUE_FORCE ; 
FUNinput.INPUTS.PLATELOADS.ISLOCAL{iSURFACE}(DOMAINS_SELECT) = 1 ; 

 