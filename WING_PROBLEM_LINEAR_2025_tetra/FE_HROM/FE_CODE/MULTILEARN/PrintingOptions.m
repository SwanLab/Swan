function DATAINM = PrintingOptions(DATAINM) ; 

DATAINM = DefaultField(DATAINM,'GID_1D_3D_print',[]); 
DATAINM.GID_1D_3D_print = DefaultField(DATAINM.GID_1D_3D_print,'ACTIVE',0); 
if DATAINM.GID_1D_3D_print.ACTIVE == 1
    GID_PRINT_3D = 0 ; 
    DATAINM.GID_1D_3D_print = DefaultField(DATAINM.GID_1D_3D_print,'METHOD_SELECTION','MAX_VON_MISES');
    DATAINM.GID_1D_3D_print = DefaultField(DATAINM.GID_1D_3D_print,'NUMBER_DOMAINS',10);
else
     GID_PRINT_3D = 1 ;  
end
DATAINM = DefaultField(DATAINM,'GID_PRINT_ALL_DOMAINS',GID_PRINT_3D) ;
DATAINM = DefaultField(DATAINM,'PRINT_GID_1D_representation',0) ;