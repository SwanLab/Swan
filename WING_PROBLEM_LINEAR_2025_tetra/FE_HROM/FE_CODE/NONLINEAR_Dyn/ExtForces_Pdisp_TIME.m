function  [FORCES,PRES_DISP,DATA] = ExtForces_Pdisp_TIME(DATA,dR,Fb,Ftrac,OPERfe,Bst)

if nargin == 0
    load('tmp.mat')
end

DATA = DefaultField(DATA,'FACTOR_TIME_BODY_FORCES',ones(size(DATA.TIME_DISCRETIZATION))) ;   % ZEROS 
DATA = DefaultField(DATA,'FACTOR_TIME_TRACTION_FORCES',DATA.TIME_DISCRETIZATION) ;  


% Body forces 
FORCES{1}.VALUE =  ConstructFequivalente(Fb,OPERfe) ; % Body forces, DOF   
FORCES{1}.VALUE_ORIGINAL = Fb ;   % Body forces  
FORCES{1}.TIME_FACTOR = DATA.FACTOR_TIME_BODY_FORCES ; %  

% Traction forces 
FORCES{2}.VALUE = ConstructFequivalente(Ftrac,OPERfe) ; ;   % Traction forces 
FORCES{2}.VALUE_ORIGINAL = Ftrac ;   % Body forces  

FORCES{2}.TIME_FACTOR = DATA.FACTOR_TIME_TRACTION_FORCES ; % 

%%% Prescribed displacements  
DATA = DefaultField(DATA,'FACTOR_TIME_DISPLACEMENTS',DATA.TIME_DISCRETIZATION) ;  
PRES_DISP.VALUE = dR ;
PRES_DISP.TIME_FACTOR = DATA.FACTOR_TIME_DISPLACEMENTS ;

% Product prescribed displacements time matrix B  (strain produced by Dir. Boundary conditions)
PRES_DISP.Bst_u = Bst(:,OPERfe.DOFs)*dR ; 


end 

function F_f  = ConstructFequivalente(Forces,OPERfe)

if ~isempty(OPERfe.DOFm)
    F_f = [Forces(OPERfe.DOFl) ; Forces(OPERfe.DOFm) + OPERfe.Gbound'*Forces(OPERfe.DOFs) ];
else
    F_f = Forces(OPERfe.DOFl)  ;
end


end

