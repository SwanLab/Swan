function GIDres2binFOLDER(INPUT_FOLDER,DELETEres,LEVEL_FOLDER,gid)


%gid = '/home/joaquin/GID_CODE_v14_01/gid' ;
disp('Converging GID ASCII res files to BINARY files ...')
% END-INPUTS
SUF = '.res' ;
SUFmsh = '.msh' ;
SUFbin = '.bin' ;
Encabezado = '#!/bin/csh -f';

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
    if ISDIR == 1 && strcmp(NAMEFILEres(1),'.')
        
    elseif ISDIR == 1 && ~strcmp(NAMEFILEres(1),'.')
        disp('***********************************************************++')
        disp(['FOLDER = ',NAMEFILEres])
        NEWINPUT_FOLDER = [INPUT_FOLDER,filesep,NAMEFILEres] ;
        % A new level of folder  is present
        [NAMETOC ]= mfilename('fullpath') ;
        OLDSCRIPT =  [NAMETOC,'.m'] ;
        NEWSCRIPT = [NAMETOC,num2str(LEVEL_FOLDER),'.m'] ;
        if exist(NEWSCRIPT,'file')
            delete(NEWSCRIPT) ;
        end
        copyfile(OLDSCRIPT,[NEWSCRIPT]) ;
        [~,NEWSCRIPTmfile,~]  = fileparts(NEWSCRIPT) ;
        LEVEL_FOLDER = LEVEL_FOLDER + 1;
        STREVAL = [NEWSCRIPTmfile,'(NEWINPUT_FOLDER,DELETEres,LEVEL_FOLDER,gid);'] ;
        eval(STREVAL) ;
        delete(NEWSCRIPT)
        
        
        
    else
        
        [~,NAMEROOT ,TERM]=  fileparts(NAMEFILEres) ;
        
        if strcmp(TERM,SUF)
            disp(['File = ',NAMEFILEres]) ;
            NAMEFILEmesh  = [NAMEROOT,SUFmsh] ;
            
            if exist([INPUT_FOLDER,filesep,NAMEFILEmesh],'file')
                
                NAMEFILEbin = [NAMEROOT,SUFbin] ;
                TXTbatch{end+1} = 'Mescape Postprocess escape' ;
                TXTbatch{end+1} = ['Mescape Files Read ','"',[INPUT_FOLDER,filesep,NAMEFILEres],'" escape'] ;
                TXTbatch{end+1} = ['Mescape Files SaveAll BinMeshesSets'] ;
                TXTbatch{end+1} = ['"',[INPUT_FOLDER,filesep,NAMEFILEbin],'" escape'] ;
                TXTbatch{end+1} = 'Mescape Quit No escape' ;
                
                FILE_BATCH = [INPUT_FOLDER,filesep,'AUX.bch'] ;
                fid = fopen(FILE_BATCH,'w') ;
                % dbstop('180')
                for i = 1:length(TXTbatch)
                    fprintf(fid,'%s\n',TXTbatch{i});
                end
                fclose(fid ) ;
                %   COMMANDFILE{end+1} = ['echo ',NAMEFILEres,' ....'] ;
                disp('Creating bin file...')
                COMMANDFILE  = [gid,' -n -b  ',FILE_BATCH] ;
                unix(COMMANDFILE) ;
                
                if DELETEres == 1
                    disp('Deleting existing files ...')
                    delete([INPUT_FOLDER,filesep,NAMEFILEres])
                    delete([INPUT_FOLDER,filesep,NAMEFILEmesh])
                end
                
            end
            
        end
    end
end

%
% FILE_BATCH = [NAMEAUX,filesep,'ALLfiles'] ;
% fid = fopen(FILE_BATCH,'w') ;
% % dbstop('180')
% for i = 1:length(COMMANDFILE)
%     fprintf(fid,'%s\n',COMMANDFILE{i});
% end
% fclose(fid ) ;
%
% disp('...Done')

%
%  unix(['chmod +x ',FILE_BATCH]) ;
%  unix([FILE_BATCH])