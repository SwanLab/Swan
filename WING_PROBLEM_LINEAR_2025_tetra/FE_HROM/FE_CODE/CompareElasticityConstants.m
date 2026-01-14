clc
clear all

FOLDER = '/home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/APUNTES_LATEX/DOCUMENTOS_anexos/MATLAB/ELASTOSTATIC_GEN/DATAWS/' ; 

NAMEDATA = { 'DATA_cuadQperio', 'DATA_cuadQmin',   'DATA_cuadQzero'}

LABELS = {'PERIODIC','MIN','ZERO','ERROR_{MIN} \%','ERROR_{ZERO} \%'} ; 
C.Ex = zeros(length(NAMEDATA),1) ; 
C.Ey = zeros(length(NAMEDATA),1) ; 
C.Ez = zeros(length(NAMEDATA),1) ; 
C.nuXY = zeros(length(NAMEDATA),1) ; 
C.nuXZ = zeros(length(NAMEDATA),1) ; 
C.nuYX = zeros(length(NAMEDATA),1) ; 
C.gYZ = zeros(length(NAMEDATA),1) ; 
C.gXZ = zeros(length(NAMEDATA),1) ; 
C.gXY = zeros(length(NAMEDATA),1) ; 



CelasGLO = cell(length(NAMEDATA),1) ; 
for i = 1:length(NAMEDATA)
load([FOLDER,'Celas_',NAMEDATA{i},'.mat']) ;
CelasGLO{i} =  Celas; 
[C.Ex(i),C.Ey(i),C.Ez(i),C.nuXY(i),C.nuXZ(i),C.nuYX(i),C.nuZY(i),C.gYZ(i),C.gXZ(i),C.gXY(i)] = ElasticityConstants(CelasGLO{i});

end

fff = fieldnames(C) ; 
VALUES = zeros(length(fff),length(NAMEDATA)) ; 
 

for i = 1:length(fff)
    for j = 1:length(NAMEDATA)
        VALUES(i,j) = C.(fff{i})(j) ; 
    end
    
     
end

VALREF = repmat(VALUES(:,1),1,length(NAMEDATA)-1) ; 
ERROR = VALUES(:,2:end) -  VALREF; 
ERROR = ERROR./VALREF*100 ; 

MATRIX = [VALUES ERROR] ; 


matrix2latex_table(MATRIX,fff,LABELS,3)

