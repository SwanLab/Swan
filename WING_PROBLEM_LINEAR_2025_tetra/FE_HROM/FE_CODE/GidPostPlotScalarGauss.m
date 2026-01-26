function GidPostPlotScalarGauss(var,nstrain,npg,nElem,time_step,NF,DATA,fid_res,ndime,setElementsRED)

% NF.VARIABLE = 'DISPLACEMENTS';
% NF.COMP(1) = 'Stress-xx';
% NF.COMP(2) = 'Stress-yy';
% NF.COMP(3) = 'Stress-zz';
% NF.COMP(4) = 'Stress-xy';
% NF.COMP(5) = 'Stress-yz';
% NF.COMP(6) = 'Stress-xz';

 if isempty(setElementsRED)
     setElements = 1:nElem ; 
 else
     setElements = setElementsRED' ;
 end
nstrain  = 1; 
LABEL_LOC = NF.VARIABLE ; %['Stress'];
% if ndime ==2
%     %         for PlainDeformationMatrix results: four components
%     % result_number_i Sxx_value Syy_value Sxy_value Szz_value
%     tipoDatAlmac = 'PlainDeformationMatrix' ;
%     nomComponente =  ['"',NF.COMP{1},'" ',NF.COMP{2},'" ',NF.COMP{3},'" ',NF.COMP{4},'" '] ; % '"StressXX" "StressYY" "StressZZ"  "StressXY" ';
% else
%     %         for Matrix results: three components ( 2D models) or six components (3D models)
%     % 2D: result_number_i Sxx_value Syy_value Sxy_value
%     % 3D: result_number_i Sxx_value Syy_value Szz_value Sxy_value Syz_value Sxz_value
%     tipoDatAlmac = 'Matrix' ;
%  %   nomComponente = '"Stress XX" "StressYY" "StressZZ" "StressXY" "StressYZ"   "StressXZ"';
%   nomComponente =  ['"',NF.COMP{1},'" ',NF.COMP{2},'" ',NF.COMP{3},'" ',NF.COMP{4},'" ',NF.COMP{5},'" ',NF.COMP{6},'" '] ;
% end
fprintf(fid_res,['Result  "',LABEL_LOC,'" "Load Analysis" ' num2str(time_step) ' Scalar OnGaussPoints "GPset"\n'],[]);
%fprintf(fid_res,'ComponentNames %s\n',nomComponente);
fprintf(fid_res,'Values\n');
format = ['%d',repmat([repmat(' %f',1,nstrain),'\n'],1,npg)];
fprintf(fid_res,format,[setElements;var']);
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