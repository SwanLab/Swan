function GidPostPlotVectNodes(d,ndime,nnod,time_step,NF,DATA,fid_res)

% NAMEFIELDS.VARIABLE = 'DISPLACEMENTS'; 
% NAMEFIELDS.COMPONENTS(1) = 'X-DISPL'; 
% NAMEFIELDS.COMPONENTS(2) = 'Y-DISPL'; 
% NAMEFIELDS.COMPONENTS(3) = 'Z-DISPL'; 


NAME = ['"',NF.VARIABLE,'"'] ; 
 


var = reshape(d,ndime,nnod) ;
%time_step = TIMEVECTOR(istep) ;
fprintf(fid_res, ['Result ',NAME,' "Load Analysis" '  num2str(time_step) ' Vector OnNodes  \n' ],[]);
if ndime==2
    fprintf(fid_res,['ComponentNames "',NF.COMP{1},'" "',NF.COMP{2},'"\n']);
elseif ndime==3
    fprintf(fid_res,['ComponentNames "',NF.COMP{1},'" "',NF.COMP{2},'" "',NF.COMP{3},'"\n']);
end
fprintf(fid_res,'Values\n');
fORMAT = ['%d',repmat(' %f',1,ndime),'\n'];
fprintf(fid_res,fORMAT,[DATA.NODES;var]);
fprintf(fid_res,'End Values\n');