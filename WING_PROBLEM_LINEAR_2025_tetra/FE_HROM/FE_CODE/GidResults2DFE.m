function GidResults2DFE(NameFile,COOR,CONNECT,TypeElement,d,strainGLO, stressGLO, ...
    React,NAME_INPUT_DATA,posgp,DATA)
% INPUT DATA
%dbstop('5')
if nargin == 0
    load('tmp2.mat')
end
%%%%%%
DATA = DefaultField(DATA,'PRINT_AVERAGE_STRESSES_ON_ELEMENTS',1) ;
elem_type = TypeElement;
NNode= size(CONNECT,2) ;;
ndime = size(COOR,2) ;
nnod = size(COOR,1) ;
npe =  size(CONNECT,2) ;
nElem = size(CONNECT,1) ;
if DATA.PRINT_AVERAGE_STRESSES_ON_ELEMENTS == 1
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
    if DATA.StrainStressWith4Components == 1
        nstrain = 4 ; 
    else
    nstrain =  3 ;
    end
else
    nstrain =6 ;
end
try
    DATA = DefaultField(DATA,'NODES',1:nnod) ;
    nnod = length(DATA.NODES) ;
catch
    DATA.NODES = 1:nnod  ;
end
stepsSH  =1; % Linear analysis, just one step
TIMEVECTOR = 1;
DATA = DefaultField(DATA,'Forces1D',[]) ;

for istepLOC = 1:length(stepsSH)
    istep = stepsSH(istepLOC) ;
    time_step = TIMEVECTOR(istep) ;
    
    %  disp(['istep=',num2str(istep)])
    % %**************************************************
    % Nodal displacements
    %**************************************************
    if ~isempty(d)
        NAMEFIELDS.VARIABLE = 'DISPLACEMENTS' ;
        NAMEFIELDS.COMP{1} = 'X-DISP' ;  NAMEFIELDS.COMP{2} = 'Y-DISP' ;  NAMEFIELDS.COMP{3} = 'Z-DISP' ;
        GidPostPlotVectNodes(d,ndime,nnod,time_step,NAMEFIELDS,DATA,fid_res) ;
    end
    
    DATA = DefaultField(DATA,'plotFLUCT',0) ;
    
    if   DATA.plotFLUCT  ==1
        d = d -DATA.dispMACRO ;
        NAMEFIELDS.VARIABLE = 'DISPL_FLUC' ;
        NAMEFIELDS.COMP{1} = 'X-FLUC' ;  NAMEFIELDS.COMP{2} = 'Y-FLUC' ;  NAMEFIELDS.COMP{3} = 'Z-FLUC' ;
        GidPostPlotVectNodes(d,ndime,nnod,time_step,NAMEFIELDS,DATA,fid_res) ;
        NAMEFIELDS.VARIABLE = 'DISPL_MACRO' ;
        NAMEFIELDS.COMP{1} = 'X-DISPmac' ;  NAMEFIELDS.COMP{2} = 'Y-DISPmac' ;  NAMEFIELDS.COMP{3} = 'Z-DISPmac' ;
        GidPostPlotVectNodes(DATA.dispMACRO,ndime,nnod,time_step,NAMEFIELDS,DATA,fid_res) ;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %**************************************************
    % Reactions
    %**************************************************
    DATA = DefaultField(DATA,'VARIABLES_PLOT_GID',[]) ; 
        DATA.VARIABLES_PLOT_GID = DefaultField(DATA.VARIABLES_PLOT_GID,'REACTIONS',1) ; 

    if ~isempty(React)  && DATA.VARIABLES_PLOT_GID.REACTIONS == 1
        NAMEFIELDS.VARIABLE = 'REACTIONS' ;
        NAMEFIELDS.COMP{1} = 'X-reac' ;  NAMEFIELDS.COMP{2} = 'Y-reac' ;  NAMEFIELDS.COMP{3} = 'Z-reac' ;
        GidPostPlotVectNodes(React,ndime,nnod,time_step,NAMEFIELDS,DATA,fid_res) ;
    end
    
    DATA = DefaultField(DATA,'Fnodes',[]);
      DATA.VARIABLES_PLOT_GID = DefaultField(DATA.VARIABLES_PLOT_GID,'Fnodes',1) ; 
      if ~isempty(DATA.Fnodes) &&   DATA.VARIABLES_PLOT_GID.Fnodes == 1
        NAMEFIELDS.VARIABLE = 'F-NODES' ;
        NAMEFIELDS.COMP{1} = 'X-force' ;  NAMEFIELDS.COMP{2} = 'Y-force' ;  NAMEFIELDS.COMP{3} = 'Z-force' ;
        GidPostPlotVectNodes(DATA.Fnodes,ndime,nnod,time_step,NAMEFIELDS,DATA,fid_res) ;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% STRESSES ........................
    % ...................................
    if ~isempty(stressGLO)
        NAMEFIELDS.VARIABLE = 'STRESS' ;
        NAMEFIELDS.COMP  = {'Stress-xx','Stress-yy','Stress-zz','Stress-xy','Stress-yz','Stress-xz'} ;
        setElementsRED = [] ;
        GidPostPlotMatGauss(stressGLO,nstrain,npg,nElem,time_step,NAMEFIELDS,DATA,fid_res,ndime,setElementsRED) ;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% REDUCED STRESSES ........................
    % ...................................
    DATA  =DefaultField(DATA,'setElementsRED',[]) ;
    if ~isempty(stressGLO) & ~isempty(DATA.setElementsRED);
        NAMEFIELDS.VARIABLE = 'GAPPY STRESS' ;
        NAMEFIELDS.COMP  = {'StressG-xx','StressG-yy','StressG-zz','StressG-xy','StressG-yz','StressG-xz'} ;
        setElementsRED = DATA.setElementsRED ;
        var = stressGLO(:,DATA.setElementsRED) ;
        GidPostPlotMatGauss(var,nstrain,npg,nElem,time_step,NAMEFIELDS,DATA,fid_res,ndime,setElementsRED) ;
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% STRAINS ........................
    % ...................................
    DATA.VARIABLES_PLOT_GID = DefaultField(DATA.VARIABLES_PLOT_GID,'STRAIN',1) ;
    if ~isempty(strainGLO)  && DATA.VARIABLES_PLOT_GID.STRAIN == 1
        NAMEFIELDS.VARIABLE = 'STRAIN' ;
        NAMEFIELDS.COMP  = {'Strain-xx','Strain-yy','Strain-zz','Strain-xy','Strain-yz','Strain-xz'} ;
        setElementsRED = [] ;
        GidPostPlotMatGauss(strainGLO,nstrain,npg,nElem,time_step,NAMEFIELDS,DATA,fid_res,ndime,setElementsRED) ;
    end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% VON MISES STRESS ........................
    % ...................................
    DATA = DefaultField(DATA,'PRINT_VONMISES_STRESS',1) ;  
    DATA = DefaultField(DATA,'stressVONMISES',[]) ;  
    if ~isempty(DATA.stressVONMISES) & DATA.PRINT_VONMISES_STRESS == 1
        NAMEFIELDS.VARIABLE = 'VON_MISES' ;       
        setElementsRED = [] ;
        GidPostPlotScalarGauss(DATA.stressVONMISES,nstrain,npg,nElem,time_step,NAMEFIELDS,DATA,fid_res,ndime,setElementsRED) ;
    end
         DATA = DefaultField(DATA,'MAXstressVONMISES',[]) ;  
    if ~isempty(DATA.MAXstressVONMISES) 
        NAMEFIELDS.VARIABLE = 'MAX_VON_MISES' ;       
        setElementsRED = [] ;
        GidPostPlotScalarGauss(DATA.MAXstressVONMISES,nstrain,npg,nElem,time_step,NAMEFIELDS,DATA,fid_res,ndime,setElementsRED) ;
    end
    
    if ~isempty(DATA.Forces1D)
        
        fff = fieldnames(DATA.Forces1D) ; 
        for ifield = 1:length(fff)
            FORCE = DATA.Forces1D.(fff{ifield}) ; 
            if ~isempty(FORCE)
                NAMEFIELDS.VARIABLE = fff{ifield} ; 
                 GidPostPlotScalarNodes(FORCE,ndime,nnod,time_step,NAMEFIELDS,DATA,fid_res) ;
            end
        end
    end
    
end




fclose(fid_res);