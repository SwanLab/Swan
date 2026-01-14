function [NameFile_msh,NameFile_res] = NameFilePostProcessGid(DATA,NAME_INPUT_DATA)


% Name of the mesh file
%if ~isempty(NameFileMesh)
%[dummy1 NameFileMeshHERE dummy2]= fileparts(NameFileMesh) ;
%else
NameFileMeshHERE = '' ;
%end

DATA = DefaultField(DATA,'POST_LabelFileGid','') ;

[DIRinp,NAMEFILEinput ]=fileparts(NAME_INPUT_DATA) ;
cdddd = cd ;
switch DIRinp
    case {'DATA_input',''}
        NameFile_msh = [cdddd,filesep,'GIDPOST',filesep,NameFileMeshHERE,NAME_INPUT_DATA,DATA.POST_LabelFileGid ,'.msh'] ;
        % Name of the results file
        NameFile_res= [cdddd,filesep,'GIDPOST',filesep,NameFileMeshHERE,NAME_INPUT_DATA,DATA.POST_LabelFileGid,'.res'] ;
        
    otherwise
        
        ROOTFOLDER = [DIRinp,filesep,'GIDPOST',filesep] ;
        if ~exist(ROOTFOLDER)
            mkdir(ROOTFOLDER) ;
        end
        NameFile_msh = [ROOTFOLDER,NameFileMeshHERE,NAMEFILEinput,DATA.POST_LabelFileGid ,'.msh'] ;
        
        NameFile_res= [ROOTFOLDER,NameFileMeshHERE,NAMEFILEinput,DATA.POST_LabelFileGid,'.res'] ;
        
        
end