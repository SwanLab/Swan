function ChangeName_mat_gid_4lev(INPUT_FOLDER)
% THIS SCRIPT SEARCH FOR ALL FILES NAMED *.gid/*.mat, and rename them to
% *./gid/*.rename_to_mat
% JAHO, 6-May-2023, Barcelona
% ----------------------------------------------------------------------------
% INPUTS
if nargin == 0
    INPUT_FOLDER = cd ;
end

% SUF_folder = '.gid' ;
% SUF_inside = '.mat' ;
%Encabezado = '#!/bin/csh -f';


% Find all files ending in SUF
TXTbatch = {} ;
AAA =dir(INPUT_FOLDER) ;

%COMMANDFILE = {} ;
%COMMANDFILE{end+1} = Encabezado ;
%NAMEAUX = [INPUT_FOLDER,filesep,'AUXFOLDER'] ;
%if exist(NAMEAUX,'dir')
%    rmdir(NAMEAUX,'s') ;
%end
%mkdir(NAMEAUX) ;

for ifile = 1:length(AAA)
    NAMEFILEres = AAA(ifile).name ;
    ISDIR = AAA(ifile).isdir ;
    TXTbatch = {} ;
    if ISDIR == 1 && length(NAMEFILEres) < 5
        % The name is too short to compare it with .gid (SUF_folder)
    else
        if ISDIR == 1 && strcmp(NAMEFILEres(end-3:end),'.gid')
            disp('***********************************************************++')
            disp(['FOLDER = ',NAMEFILEres])
            % Rename .mat
            NameFolderLOC = [INPUT_FOLDER,'/',NAMEFILEres] ;
            AAAfiles =dir(NameFolderLOC) ;
            for ifile = 1:length(AAAfiles)
                NAMEFILEgid = AAAfiles(ifile).name ;
                if  length(NAMEFILEgid) < 5
                else
                    if  strcmp(NAMEFILEgid(end-3:end),'.mat')
                        NameOld = [NameFolderLOC,'/',NAMEFILEgid(1:end-4),'.mat'] ; ;
                        NameNew = [NameFolderLOC,'/',NAMEFILEgid(1:end-4),'.rename_2_mat'] ; ;
                        disp([NAMEFILEgid(1:end-4),'.mat','--->',NAMEFILEgid(1:end-4),'.rename_2_mat']);
                        unix([ ' cp ',NameOld,' ',NameNew]);
                    end
                end
                %
                %              if ~strcmp(old_name(end),'/') & ~strcmp(new,old)
                %         disp([old,'--->',new]);
                %         unix(['cd ',PathProj,' ; mv ',old,' ',new]);
            end
        elseif ISDIR == 1  
            NameFolderLOC = [INPUT_FOLDER,'/',NAMEFILEres] ;
            ChangeName_mat_gid_3lev(NameFolderLOC)  ; 
        end
    end
    
    
    
    
end
