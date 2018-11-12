function exportSTL (fname,algorithm,user_paths)

%% Initialize
down_path = user_paths.temp_folder;
[stlname,path] = uiputfile('*.stl','Save STL as',[user_paths.default_stl_path,'/']);
if ~path
    return;
end

% Copy file to temporal path
[org_path,file,ext] = fileparts(fname);
res_file = [down_path,'/',file,ext];
msh_file = [down_path,'/',file,'.msh'];
copyfile(fname,down_path,'f');
fname = [org_path,'/',file,'.msh'];
copyfile(fname,down_path,'f');

% Function to fix directions
fixdir = @(str) strrep(str,'\','/');

% Auxiliar files
ascii_file = [user_paths.temp_folder,'/ascii'];
geom_file = [user_paths.temp_folder,'/geom.dxf'];

%% Run GID
if strcmp(algorithm,'level_set')
    macro_name = 'GeomSTL';
else
    macro_name = 'GeomSTL_PG';
end

fix_macros_GID(user_paths);
if ispc % WINDOWS
    command = sprintf('start /wait %s -t "::toolbarmacros::macrospace::%s  "%s" "%s" "%s" "', ...
                        user_paths.gid_path,macro_name,fixdir(res_file),fixdir(ascii_file),fixdir(geom_file));
    dos(command);
elseif isunix % LINUX
    command = sprintf('%s -t "::toolbarmacros::macrospace::%s  "%s" "%s" "%s" "', ...
                        user_paths.gid_path,macro_name,fixdir(res_file),fixdir(ascii_file),fixdir(geom_file));
    unix(command);

else
    error('Operating system not included!.');
end

%% Delete temporal files
delete(res_file);
delete(msh_file);
files = ls(down_path);
files = cellstr(files);
files(1:2) = [];

del = {'ascii' 'geom'};
N = length(files);
for i = 1:N
    [~,fname,~] = fileparts(files{i});
    if ismember(fname,del)
        delete([down_path,'/',files{i}]);
    end
end

%% Export final STL mesh
stl = [down_path,'/result.stl'];
movefile(stl,[path,stlname],'f')

end

