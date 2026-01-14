function GidResults2DFEDyn(NameFile,COOR,CONNECT,TypeElement,d,strainGLO, stressGLO,  React,NAME_INPUT_DATA,posgp,t)
% INPUT DATA
%dbstop('4')
if nargin == 0
    load('tmp1.mat')
end
%%%%%%
elem_type = TypeElement;
NNode= size(CONNECT,2) ;;
ndime = size(COOR,2) ;
nnod = size(COOR,1) ;
npe =  size(CONNECT,2) ;
nElem = size(CONNECT,1) ;
npg = size(posgp,2);
fid_res = fopen(NameFile,'wt');
fprintf(fid_res,'GiD Post Results File 1.0 \n');
xg =posgp ;
fprintf(fid_res,['GaussPoints "GPset" Elemtype ',elem_type,'\n']);
fprintf(fid_res,['Number of Gauss Points: ',num2str(npg),'\n']);
fprintf(fid_res,'Nodes not included\n');
fprintf(fid_res,'Natural Coordinates: Internal\n');
fprintf(fid_res,'End GaussPoints\n');
if ndime == 2
    nstrain =4 ;
else
    nstrain =6 ;
end
stepsSH  =1:length(t); %  
TIMEVECTOR = t;
for istepLOC = 1:length(stepsSH)
    istep = stepsSH(istepLOC) ;
    disp(['istep=',num2str(istep)])
    % %**************************************************
    % Nodal displacements
    %**************************************************
    var = reshape(d(:,istepLOC),ndime,nnod) ;
    time_step = TIMEVECTOR(istep) ;
    fprintf(fid_res, ['Result "Nodal displacement"  "Load Analysis" '  num2str(time_step) ' Vector OnNodes  \n' ],[]);
    if ndime==2
        fprintf(fid_res,'ComponentNames "X-DISPL" "Y-DISPL"\n');
    elseif ndime==3
        fprintf(fid_res,'ComponentNames "X-DISPL" "Y-DISPL" "Z-DISPL"\n');
    end
    fprintf(fid_res,'Values\n');
    fORMAT = ['%d',repmat(' %f',1,ndime),'\n'];
    fprintf(fid_res,fORMAT,[1:nnod;var]);
    fprintf(fid_res,'End Values\n');
%     if ~isempty(React)
%         var = reshape(React,ndime,nnod) ;
%         time_step = TIMEVECTOR(istep) ;
%         fprintf(fid_res, ['Result "Nodal reactions"  "Load Analysis" '  num2str(time_step) ' Vector OnNodes  \n' ],[]);
%         if ndime==2
%             fprintf(fid_res,'ComponentNames "X-REAC" "Y-REAC"\n');
%         elseif ndime==3
%             fprintf(fid_res,'ComponentNames "X-REAC" "Y-REAC" "Z-REAC"\n');
%         end
%         fprintf(fid_res,'Values\n');
%         fORMAT = ['%d',repmat(' %f',1,ndime),'\n'];
%         fprintf(fid_res,fORMAT,[1:nnod;var]);
%         fprintf(fid_res,'End Values\n');
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%% STRESSES ........................
%     % ...................................
%     var = stressGLO ;
%     LABEL_LOC = ['Stress'];
%     if ndime ==2
%         %         for PlainDeformationMatrix results: four components
%         % result_number_i Sxx_value Syy_value Sxy_value Szz_value
%         tipoDatAlmac = 'PlainDeformationMatrix' ;
%         nomComponente = '"StressXX" "StressYY" "StressZZ"  "StressXY" ';
%     else
%         %         for Matrix results: three components ( 2D models) or six components (3D models)
%         % 2D: result_number_i Sxx_value Syy_value Sxy_value
%         % 3D: result_number_i Sxx_value Syy_value Szz_value Sxy_value Syz_value Sxz_value
%         tipoDatAlmac = 'Matrix' ;
%         nomComponente = '"Stress XX" "StressYY" "StressZZ" "StressXY" "StressYZ"   "StressXZ"';
%     end
%     fprintf(fid_res,['Result  "',LABEL_LOC,'" "Load Analysis" ' num2str(time_step) ' Matrix OnGaussPoints "GPset"\n'],[],tipoDatAlmac);
%     fprintf(fid_res,'ComponentNames %s\n',nomComponente);
%     fprintf(fid_res,'Values\n');
%     format = ['%d',repmat([repmat(' %f',1,nstrain),'\n'],1,npg)];
%     
%     %Se est� sacando el �ltimo valor del vector de variable de hist�rica de cada punto de gauss, que es
%     %justamente la deformaci�n pl�stica equivalente.
%     fprintf(fid_res,format,[1:nElem;var]);
%     fprintf(fid_res,'End Values\n');
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%% STRAINS ........................
%     % ...................................
%     var = strainGLO ;
%     LABEL_LOC = ['Strain'];
%     if ndime ==2
%         %         for PlainDeformationMatrix results: four components
%         % result_number_i Sxx_value Syy_value Sxy_value Szz_value
%         tipoDatAlmac = 'PlainDeformationMatrix' ;
%         nomComponente = '"StrainXX" "StrainYY" "StrainZZ"  "StrainXY" ';
%     else
%         %         for Matrix results: three components ( 2D models) or six components (3D models)
%         % 2D: result_number_i Sxx_value Syy_value Sxy_value
%         % 3D: result_number_i Sxx_value Syy_value Szz_value Sxy_value Syz_value Sxz_value
%         tipoDatAlmac = 'Matrix' ;
%         nomComponente = '"Strain XX" "StrainYY" "StrainZZ" "StrainXY" "StrainYZ"   "StrainXZ"';
%     end
%     fprintf(fid_res,['Result  "',LABEL_LOC,'" "Load Analysis" ' num2str(time_step) ' Matrix OnGaussPoints "GPset"\n'],[],tipoDatAlmac);
%     fprintf(fid_res,'ComponentNames %s\n',nomComponente);
%     fprintf(fid_res,'Values\n');
%     format = ['%d',repmat([repmat(' %f',1,nstrain),'\n'],1,npg)];
%     
%     %Se est� sacando el �ltimo valor del vector de variable de hist�rica de cada punto de gauss, que es
%     %justamente la deformaci�n pl�stica equivalente.
%     fprintf(fid_res,format,[1:nElem;var]);
%     fprintf(fid_res,'End Values\n');
    
    
    
end



fclose(fid_res);