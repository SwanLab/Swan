function GidResults2DFE_modesGAUSS(NameFile,COOR,CONNECT,TypeElement,MODES,posgp,DATA)
% INPUT DATA
%dbstop('4')
if nargin == 0
    load('tmp.mat')
end


COMP{1} = 'Stress-xx';
COMP{2} = 'Stress-yy';
COMP{3} = 'Stress-zz';
COMP{4} = 'Stress-xy';
COMP{5} = 'Stress-yz';
COMP{6} = 'Stress-xz';


%%%%%%
elem_type = TypeElement;
NNode= size(CONNECT,2) ;;
ndime = size(COOR,2) ;
nnod = size(COOR,1) ;
npe =  size(CONNECT,2) ;
nElem = size(CONNECT,1) ;
setElements = 1:nElem ;
if isempty(posgp)
    npg = 1;
else
    npg = size(posgp,2);
end
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



stepsSH  =1; % Linear analysis, just one step
TIMEVECTOR = 1;
istep = 1;
disp(['istep=',num2str(istep)])
% %**************************************************
% Natural modes
%**************************************************
%DATA.NODES = NODESref ;
%dbstop('40')
% try
%     DATA = DefaultField(DATA,'NODES',1:nnod) ;
%     nnod = length(DATA.NODES) ;
% catch
%     DATA.NODES = 1:nnod  ;
% end

if ndime ==2
    %         for PlainDeformationMatrix results: four components
    % result_number_i Sxx_value Syy_value Sxy_value Szz_value
    tipoDatAlmac = 'PlainDeformationMatrix' ;
    nomComponente =  ['"',COMP{1},'" "',COMP{2},'" "',COMP{4},'" "',COMP{3},'" '] ; % '"StressXX" "StressYY" "StressZZ"  "StressXY" ';
else
    %         for Matrix results: three components ( 2D models) or six components (3D models)
    % 2D: result_number_i Sxx_value Syy_value Sxy_value
    % 3D: result_number_i Sxx_value Syy_value Szz_value Sxy_value Syz_value Sxz_value
    tipoDatAlmac = 'Matrix' ;
    %   nomComponente = '"Stress XX" "StressYY" "StressZZ" "StressXY" "StressYZ"   "StressXZ"';
    nomComponente =  ['"',COMP{1},'" "',COMP{2},'" "',COMP{3},'" "',COMP{4},'" "',COMP{5},'" "',COMP{6},'" '] ;
end


for imode = 1:size(MODES,2)
    d = MODES(:,imode) ;
    
    
    ncomp =  nstrain*size(xg,2);
    var = reshape(d,ncomp,[]) ;
    
    
    time_step = TIMEVECTOR(istep) ;
    %  fprintf(fid_res, ['Result "Mode ',num2str(imode),'"  "Load Analysis" '  num2str(time_step) ' Vector OnNodes  \n' ],[]);
    
    fprintf(fid_res,['Result "Mode',num2str(imode),'"  "Load Analysis" ', num2str(time_step),' ',tipoDatAlmac,' OnGaussPoints "GPset"\n'],[],tipoDatAlmac);
    fprintf(fid_res,'ComponentNames %s\n',nomComponente);
    fprintf(fid_res,'Values\n');
    
    format = ['%d',repmat([repmat(' %f',1,nstrain),'\n'],1,npg)];
    fprintf(fid_res,format,[setElements;var]);
    fprintf(fid_res,'End Values\n');
    %
    %     fprintf(fid_res,'Values\n');
    %     fORMAT = ['%d',repmat(' %f',1,ndime),'\n'];
    %     fprintf(fid_res,fORMAT,[DATA.NODES;var]);
    %     fprintf(fid_res,'End Values\n');
    
    
    
    
    
    
    
    
end


fclose(fid_res);