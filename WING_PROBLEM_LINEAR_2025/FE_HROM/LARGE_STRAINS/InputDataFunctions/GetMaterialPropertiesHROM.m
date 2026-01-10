function [MATPRO,INICOND] = GetMaterialPropertiesHROM(ECMdata,DATAFE,MATPRO,DATAHROM,OTHER_output) 

if nargin == 0
    load('tmp3.mat')
end

ngausT = DATAHROM.MESH.ngausT ; 
nelem =  DATAHROM.MESH.nelem ; 
nTOTALdofsGAUSS = DATAHROM.MESH.ngausT*DATAHROM.MESH.nstrain ; 



setPointsElement = [] ; 
ECMdata = DefaultField(ECMdata,'setPoints',[]) ; 

if  isempty(ECMdata.setPoints)
    % The continuous ECM has been used 
    % We assume that, for a given element, the elastic properties are the
    % same for all the Gauss points of such elements 
    setPointsElement = ECMdata.setElements*DATAFE.MESH.ngaus_STRESS; % Last Gauss point  of an element
    % For instance, if ngaus_STRESS = 9, and setElement = 1, then 
    % setPointsElement=9 
    setIndices = small2large(setPointsElement,DATAFE.MESH.nstrain) ;
   % MATPRO.celasglo = MATPRO.celasglo(setIndicesLOC,:) ;
    
else
    setIndices = small2large(ECMdata.setPoints,DATAFE.MESH.nstrain) ;
  %  MATPRO.celasglo = MATPRO.celasglo(setIndices,:) ;
    
end

% MAterial properties 
fff = fieldnames(MATPRO); 
MATPRO = GetVariablesAtECMpoints(MATPRO,fff,ngausT,ECMdata,setPointsElement,nTOTALdofsGAUSS,setIndices) ; 


 
% Initial conditions list internal variables 
 DATAFE = DefaultField(DATAFE,'ListFieldInternalVariables',[] ); 
 OTHER_output = DefaultField(OTHER_output,'INICOND',[]) ; 
 INICOND =  OTHER_output.INICOND; 
 if ~isempty(DATAFE.ListFieldInternalVariables)
    INICOND = GetVariablesAtECMpoints(INICOND,DATAFE.ListFieldInternalVariables,ngausT,ECMdata,setPointsElement,nTOTALdofsGAUSS,setIndices) ;
 
     
 end

 