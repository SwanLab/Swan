function CreateBatchFileGid_mesh_data(FILE,GID_FOLDER)

if nargin == 1
    GID_FOLDER = '/home/joaquin/GID_CODE_v13/' ;
end

% disp('Mescape Meshing Generate')
% disp('Yes')
% disp('2 MeshingParametersFrom=Preferences')
% disp('Mescape Meshing MeshView')
% disp('Mescape Meshing CreateBoundary')
% disp('Yes')
% disp(['Mescape'])



disp(['Mescape Files WriteMesh'])
disp([FILE,'.msh'])
disp(['Files WriteCalcFile'])
disp([FILE,'.dat'])
disp(['Mescape'])

WRITE_BATCH = 0; 

if WRITE_BATCH ==1 

TXT = {} ;
TXT{end+1} = 'Mescape Meshing MeshView';
TXT{end+1} = 'Mescape Meshing CreateBoundary';
TXT{end+1} = 'Yes';
TXT{end+1} = 'escape';

TXT{end+1} = ['Mescape Files WriteMesh'];
TXT{end+1} = [FILE,'.msh'];
TXT{end+1} = ['Files WriteCalcFile'];
TXT{end+1} = [FILE,'.dat'];
TXT{end+1} =  ['Mescape'];
 TXT{end+1} = 'Mescape Quit No escape' ;


[DIRL,~,~] = fileparts(FILE) ;

FILE_BATCH = [DIRL,filesep,'AUX.bch'] ;
fid = fopen(FILE_BATCH,'w') ;
% dbstop('180')
for i = 1:length(TXT)
    fprintf(fid,'%s\n',TXT{i});
end
fclose(fid ) ;

disp('Automatically creating boundary mesh and DATs GID  files...')
COMMANDFILE  = [GID_FOLDER,filesep,'gid -n -b  ',FILE_BATCH] ;
unix(COMMANDFILE) ;

end


