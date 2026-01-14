function DATAOUT = mainELASTOSTATIC(FUNinput)
format long g
% Finite Element Program for Elastostatic problems  
% ECA.
% Technical University of Catalonia
% JoaquIn A. Hdez, October 23-th, 2015
% ---------------------------------------------------
if exist('ElemBnd')==0
    addpath('ROUTINES_AUX') ;
end

 
%%% INPUT  %%% 
% PREPROCESS  
[COOR,CN,TypeElement,TypeElementB, celasglo,  DOFr,dR,...  
    Tnod,CNb,fNOD,Fpnt,NameFileMesh,typePROBLEM,celasgloINV,CONNECTb,DOFm,Gb,...
    DATA,MaterialType] = ReadInputDataFile(FUNinput)  ; 

% SOLVER 
% --------------------------------------------
 [d strainGLO stressGLO  React posgp CN MaterialType]= SolveElastFE(COOR,CN,TypeElement,TypeElementB,...
    celasglo,  DOFr,dR,Tnod,CNb,fNOD,Fpnt,typePROBLEM,celasgloINV,DATA,CONNECTb,DOFm,Gb,MaterialType)  ; 

% POSTPROCESS
% --------------------------------------------
GidPostProcess(COOR,CN,TypeElement,d,strainGLO, stressGLO,  React,NAME_INPUT_DATA,posgp,NameFileMesh,MaterialType,DATA);