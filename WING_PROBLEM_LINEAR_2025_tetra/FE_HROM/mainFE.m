function DATAOUT = mainFE(INPUTDATA,DATA)
% Finite element code for composite materials
%  Technical University of Catalonia, DIC-2015/NOV-2016/APRIL-2018
%  % Joaquín A. HERNANDEZ, jhortega@cimne.upc.edu
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-
if nargin == 0
    INPUTDATA =  'Dplate_test2' ; 'DATA_BEAM3D' ;   'DplateFINEMESH' ;'DATA_THIN_XY' ; 'DATA_POROUS_left' ; 'DATA_BEAM3DrepSdomDF30' ;  'DATA_BEAM3DrepSdomPt' ; 'DATA_BEAM3DrepSdomPz' ;  'DATA_BEAM3DrepTEST6' ; 'DATA_BEAM3Drep1domPP' ;   'DATA_BEAM3DrepMM2' ; 'DATA_BEAM3Drep' ; 'DATA_BEAM3DrepMM2' ;'DATA_BEAM3DrepMM' ;'DATA_BEAM3Drep' ;'DATA_BEAM3Drep' ;   'DATAStent1';'DATAStent1';'DATA_BEAMhexaINDt3';   'DATA_BEAMthinw';   'DATA_BEAMsolid' ; 'DATA_BEAMhexat';  'DATA_BEAMporousREPt' ;  'DATA_BEAMslice' ;    'DATA_BEAM_sten'; 'DATA_BEAMporousBEND'; 'DATA_BEAM_Hexa2D_D6_bend';'DATA_BEAM3dD_aeroCOMP' ; 'DATA_BEAM_Hexa2D_D6_bend'; 'DATA_BEAMporousBENDaxial' ; 'DATA_BEAM_Hexa2D_bend' ;  'DATA_BEAM_Hexa2D_gen' ;  'DATA_BHexaD4_INDENT' ;  'DATA_BHexa2_invINDENT'; 'DATA_HexaD4compr' ;  'DATA_BHexa2_COMPR' ;  'DATA_BHexa2hom_COMPR' ;   'DATA_BHexa2hom_invINDENT';  'DATA_BHexa2hom_INDENT' ;'DATA_BHexa2_INDENT' ;   'DATA_BHexa2_invINDENT' ;'DATA_BHexa2hom_invINDENT'  ; 'DATA_BHexa2_INDENT2' ; 'DATA_BEAM_Hexa2D_gen' ;'DATA_BEAM_Hexa2D_HOMOG_GEN'  ;   ;'DATA_BHexa2_INDENT' ; 'DATA_BHexa2hom_COMPR';
end
%addpath('DATA_input')
if ~exist('Quadratic3N')
    addpath(genpath(getenv('MATLAB_CODES_FEHROM')))  ;
end
%DATA.PLAY_SOUND = 1; 

% Check if the folder in which INPUTDATA is defined is on the current
% folder

[DIRinp,NAMEFILEinput ]=fileparts(INPUTDATA) ; 

addpath(DIRinp) ; 
eval(NAMEFILEinput) ;
DATA.INPUTDATAfile = INPUTDATA ;
DATA.typePROBLEM  = typePROBLEM ; 
rmpath(DIRinp) ; 

DATA = DefaultField(DATA,'NameFileMeshDATA',FUNinput.INPUTS.NameFileMesh) ; 
FUNinput.INPUTS.NameFileMesh = DATA.NameFileMeshDATA ; 

%----------------------------------------------------
% Calling Finite Element elastostatic program
% ---------------------------------------------------
%----------------------------------------------------
tic
DATAOUT = FE_ELASTOSTATIC(FUNinput,DATA) ;
toc
% ----------------------------------------------------
DATA = DefaultField(DATA,'PLAY_SOUND',0) ; 
if DATA.PLAY_SOUND ==1
load handel.mat;
sound(y, Fs);
end


