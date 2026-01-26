function [DATAOUT DATA]= FE_ELASTOSTATIC(FUNinput,DATA)
format long g
% Finite Element Program for Elastostatic problems
% Technical University of Catalonia
% JoaquIn A. Hdez, October 23-th, 2015
% ----------------------------------------------------------
% INPUTS:
% ------------------------------------
%%% Structure array: FUNinput
% ----------------------------
% FUNinput.NAME  --> Name of the function employed to automatically
% generate geometric, material, and boundary conditions data
% Available options: 1. "INPUT_PERIODIC"
% FUNinput.INPUTS : Data structure containing the input data for function
% " FUNinput.NAME"
% For instance, for FUNinput.NAME = "INPUT_PERIODIC", the input data are
%%%%  FUNinput.INPUTS.strainINP = Macroscopic strain
%%%%  FUNinput.INPUTS.MATERIAL =  MAterial properties
%%%%  FUNinput.INPUTS.NameFileMesh =  Name of the file with mesh data
%%%% Structure array: DATA
% Miscellaneous data
%
% ---------------------------------------------------
% if exist('ElemBnd')==0
%     addpath('ROUTINES_AUX') ;
% end
if exist('INPUT_PERIODIC')==0
    addpath('FE_CODE') ;
end
if exist('GeometryStructure','file') == 0 ; addpath('FE_CODE/BeamROM') ; end


%%% INPUT  %%%
% PREPROCESS
%dbstop('32')
[COOR,CN,TypeElement,TypeElementB, celasglo,  DOFr,dR,...
    Tnod,CNb,fNOD,Fpnt,NameFileMesh,typePROBLEM,celasgloINV,CONNECTb,DOFm,Gb,...
    DATA,MaterialType,DOMAINVAR] = ReadInputDataFile(FUNinput,DATA)  ;

% SOLVER
% --------------------------------------------
DATA  = DefaultField(DATA,'NOCALCULATE_DISPLACEMENTS',0) ;
[d, strainGLOgid, stressGLOgid,    React, posgp, CN, MaterialType, DATAOUT,Fnodes]=...
    SolveElastFE(COOR,CN,TypeElement,TypeElementB,...
    celasglo,  DOFr,dR,Tnod,CNb,fNOD,Fpnt,typePROBLEM,celasgloINV,DATA,CONNECTb,DOFm,Gb,MaterialType)  ;




DATAOUT.DOMAINVAR = DOMAINVAR ;
DATAOUT.disp = d ;
DATAOUT.posgp = posgp ;
DATAOUT.React = React ;
DATAOUT.Fnodes = Fnodes ;

DATAOUT.nelem = size(CN,1) ;
try
    DATA.stressVONMISES = DATAOUT.stressVONMISES;
catch
    DATA.stressVONMISES = [] ;
end

% POSTPROCESS
% --------------------------------------------
if isfield(DATA,'INPUTDATAfile')
    NAME_BASE_GIDfiles = DATA.INPUTDATAfile ;
    
else
    
    NAME_BASE_GIDfiles = FUNinput.NAME ;
end
if  DATA.NOCALCULATE_DISPLACEMENTS == 0
    %     DATA = DefaultField(DATA,'INPUTDATAfile',[]) ;
    %     if ~isempty(DATA.INPUTDATAfile )
    %         [DIRinp,NAMEFILEinput ]=fileparts(DATA.INPUTDATAfile) ;
    %         addpath(DIRinp) ;
    %         eval(NAMEFILEinput ) ;
    %         rmpath(DIRinp) ;
    %
    %
    %         DATA.MATERIAL = MATERIAL ;
    %     else
    %         DATA.MATERIAL =[]  ;
    %     end
    if isfield(DATA,'MATERIAL_ORIGINAL')
        DATA.MATERIAL = DATA.MATERIAL_ORIGINAL ;
    end
    % DATA.MATERIAL_ORIGINAL = INPUTS_LOC.MATERIAL ;
    DATA = DefaultField(DATA,'PRINT_GID',[]) ;
    DATA.PRINT_GID = DefaultField(DATA.PRINT_GID,'PRINT_RESULTS',1) ;
    if DATA.PRINT_GID.PRINT_RESULTS == 1
        GidPostProcess(COOR,CN,TypeElement,d,strainGLOgid, stressGLOgid,  ...
            React,NAME_BASE_GIDfiles,posgp,NameFileMesh,MaterialType,DATA,Fnodes);
    else
        DATAOUT.React =  sparse(DATAOUT.React)  ;
        DATAOUT.Fnodes = sparse(DATAOUT.Fnodes) ;
        disp('...Gid post-process file was not printed')
    end
% else
%     DATA = DefaultField(DATA,'NOCALCULATE_DISPLACEMENTS_BUT_PRINT',0) ; 
%     if DATA.NOCALCULATE_DISPLACEMENTS_BUT_PRINT ==1
%         
%         strainGLOgid = [] ; stressGLOgid = [] ; React = [] ; Fnodes = [] ; d = [] ; 
%         
%         GidPostProcess(COOR,CN,TypeElement,d,strainGLOgid, stressGLOgid,  ...
%             React,NAME_BASE_GIDfiles,posgp,NameFileMesh,MaterialType,DATA,Fnodes);
%     end
end