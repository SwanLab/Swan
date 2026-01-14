
disp('*******************************************************+')
disp('Freeing memory ...')
if DATA.MINIMIZE_MEMORY_REQUIREMENTS_IN_SOLVING == 1
LIST_VARIABLES_TO_DELETE = {'celasglo','celasgloINV','Nst','Cglo','Bst'};
elseif DATA.MINIMIZE_MEMORY_REQUIREMENTS_IN_SOLVING == 3
    LIST_VARIABLES_TO_DELETE = {'celasglo','celasgloINV','Nst'};
else
    error('Option not implemented')
end
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


Nst = [] ; celasglo = [] ; celasgloINV = []; 