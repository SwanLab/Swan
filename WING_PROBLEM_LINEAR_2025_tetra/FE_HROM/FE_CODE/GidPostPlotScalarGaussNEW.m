function GidPostPlotScalarGaussNEW(VAR,npg,setElements,time_step,NF,fid_res,NAMEMESH) ;
 
%VAR = reshape(VAR,nstrain,length(setElements)) ;

nstrain  = 1; 
VAR = reshape(VAR,npg,length(setElements)) ;
LABEL_LOC = NF.VARIABLE ;
fprintf(fid_res,['Result  "',LABEL_LOC,'" "Load Analysis" ' num2str(time_step) ' Scalar OnGaussPoints "',NAMEMESH,'"\n'],[]);
%fprintf(fid_res,'ComponentNames %s\n',nomComponente);
fprintf(fid_res,'Values\n');
format = ['%d',repmat([repmat(' %f',1,nstrain),'\n'],1,npg)];
fprintf(fid_res,format,[setElements;VAR]);
fprintf(fid_res,'End Values\n');



%
%
% NAME = ['"',NAMEFIELDS.VARIABLE,'"'] ;
%
%
% var = reshape(d,ndime,nnod) ;
% time_step = TIMEVECTOR(istep) ;
% fprintf(fid_res, ['Result ',NAME,' "Load Analysis" '  num2str(time_step) ' Vector OnNodes  \n' ],[]);
% if ndime==2
%     fprintf(fid_res,'ComponentNames "X-DISPL" "Y-DISPL"\n');
% elseif ndime==3
%     fprintf(fid_res,'ComponentNames "X-DISPL" "Y-DISPL" "Z-DISPL"\n');
% end
% fprintf(fid_res,'Values\n');
% fORMAT = ['%d',repmat(' %f',1,ndime),'\n'];
% fprintf(fid_res,fORMAT,[DATA.NODES;var]);
% fprintf(fid_res,'End Values\n');