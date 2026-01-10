function GidResults2DFE_modesFINT(NameFile,COOR,CONNECT,TypeElement,MODES,posgp,DATA)
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
% DATA = DefaultField(DATA,'NODES',1:nnod) ; 
% nnod = length(DATA.NODES) ; 
% catch
%     DATA.NODES = 1:nnod  ; 
% end
setElements = 1:nElem ; 
NAMEMESH = 'GPset' ; 


DATA= DefaultField(DATA,'LABEL_NAME_VARIABLE','Fint') ; 
LABEL = DATA.LABEL_NAME_VARIABLE ; 
 

for imode = 1:size(MODES,2)
  
    VAR = MODES(:,imode); 
    nstrain = 1; 
    time_step = 1; 
    NF.VARIABLE = ['Mode ',LABEL,'=',num2str(imode)] ; 
    
    GidPostPlotScalarGaussNEW(VAR,npg,setElements,time_step,NF,fid_res,NAMEMESH) ;
    
    
     
    
    
end


fclose(fid_res);