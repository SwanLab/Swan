
disp('*******************************************************+')
disp('Freeing memory ...')
LIST_VARIABLES_TO_DELETE = {'celasglo','celasgloINV','Nst','Cglo','Bst'};
[NWS] = fileparts(DATA.nameWORKSPACE);
NWS_Cglo = [NWS,filesep,'Cglo.mat'] ; 
NWS_Bst = [NWS,filesep,'Bst.mat'] ; 
disp('Saving Cglo...')
save(NWS_Cglo,'Cglo','-v7.3') ;
disp('Saving Bst...')
save(NWS_Bst,'Bst','-v7.3') ;

DATA.nstrain =  size(celasglo,1) ;


[AAAall] = whos ;
[AAA] = whos(LIST_VARIABLES_TO_DELETE{:});
bytesDELETE = 0 ;
for iii = 1:length(AAA)
    bytesDELETE = bytesDELETE + AAA(iii).bytes ;
end
bytesDELETE  =bytesDELETE*1e-9;
bytesALL = 0 ; 
for iii = 1:length(AAAall)
    bytesALL = bytesALL + AAAall(iii).bytes ;
end
bytesALL = bytesALL*1e-9 ;
clear(LIST_VARIABLES_TO_DELETE{:});
peprc = bytesDELETE/bytesALL*100; 
disp(['Deleted =',num2str(peprc),' %  (',num2str(bytesDELETE),' Gb  of ',  num2str(bytesALL),' Gb)'])
disp('*******************************************************+')
