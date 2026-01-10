function [DATAOUT, DATA,OPERfe,NODES_SNAP]= FE_NONLINEAR(FUNinput,DATA)
format long g
% Finite Element Program for Nonlinear Dynamic problems. Copy of
% FE_ELASTOSTATIC.m
% Technical University of Catalonia
% JoaquIn A. Hdez, September 18-th, 2018
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
if nargin == 0
    load('tmp.mat')
end

if exist('INPUT_PERIODIC')==0
    addpath('FE_CODE') ;
end

if exist('GeometryStructure','file') == 0 ; addpath('FE_CODE/BeamROM') ; end


%%% INPUT  %%%
% PREPROCESS
%dbstop('32')
[COOR,CN,TypeElement,TypeElementB, PROPMAT,  DOFr,dR,...
    Tnod,CNb,fNOD,Fpnt,NameFileMesh,typePROBLEM,~,CONNECTb,DOFm,Gb,...
    DATA,MaterialType,DOMAINVAR] = ReadInputDataFile(FUNinput,DATA)  ;

% SOLVER
% --------------------------------------------
DATA  = DefaultField(DATA,'NOCALCULATE_DISPLACEMENTS',0) ;
[NODES_SNAP,GAUSS_SNAP,OPERfe,DATA,DATAOUT,NODESV_PROP,GAUSSV_PROP]=...
    SolveNonLinearFE(COOR,CN,TypeElement,TypeElementB,...
    PROPMAT,  DOFr,dR,Tnod,CNb,fNOD,Fpnt,typePROBLEM,DATA,CONNECTb,DOFm,Gb,MaterialType)  ;
% DATAOUT.DOMAINVAR = DOMAINVAR ;
% DATAOUT.disp = d ;
% DATAOUT.posgp = posgp ;
% DATAOUT.React = React ;
% try
%     DATA.stressVONMISES = DATAOUT.stressVONMISES;
% catch
%     DATA.stressVONMISES = [] ;
% end

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
    
    disp('Printing results in GID')
    DATA.PROPMAT = PROPMAT;
    DATA = DefaultField(DATA,'PRINT_AVERAGE_GAUSS_VARIABLES',0) ;  % Pending to be implemented

    GidPostProcess_GEN(COOR,CN,TypeElement,DATA,NAME_BASE_GIDfiles,NameFileMesh,MaterialType,...
        NODES_SNAP,GAUSS_SNAP,OPERfe,NODESV_PROP,GAUSSV_PROP);
    disp('...End')
    
    
    %
  
    
  
    
end