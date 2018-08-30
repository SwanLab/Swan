function fix_macros_GID (user_paths)

%% Read TCL into cell
fixdir = @(str) strrep(str,'\','/');
fname = user_paths.macro_dest;
fid = fopen(fname,'r');

i = 1;
tline = fgetl(fid);
data{i,1} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    data{i,1} = tline;
end
fclose(fid);

%% Export STL functions
% Arguments modification
pattern = 'proc Geom';
k = regexpi(data(1:end-1),pattern);
index=~cellfun('isempty', k);

args = '{res_file ascii_file geom_file}';
data(index) = add_string(data(index),args,'{}');

% DXF path modification
pattern = 'GiD_Process Mescape Files WriteForBAS';
k = regexpi(data(1:end-1),pattern);
index= ~cellfun('isempty', k);
index = filter_rang(data,index,'proc Geom');
data(index) = {['    GiD_Process Mescape Files WriteForBAS {',user_paths.dxf_path,'/DXF.bas} ${geom_file}']};

% STL path modification
pattern = 'GiD_Process Mescape Files WriteForBAS';
k = regexpi(data(1:end-1),pattern);
index= ~cellfun('isempty', k);
index = filter_rang(data,index,'proc Export');
data(index) = {['    GiD_Process Mescape Files WriteForBAS {',user_paths.dxf_path,'/STL.bas} {',fixdir(user_paths.temp_folder),'/result.stl}']};

%% Export PNG functions
pattern = 'proc PNG';
k = regexpi(data(1:end-1),pattern);
index=~cellfun('isempty', k);

args = '{res_file view_file png_file corner1 corner2 corner3 corner4 corner5 corner6 symcase}';
data(index) = add_string(data(index),args,'{}');

%% Write cell back into file
fid = fopen(fname, 'wt');
for i = 1:numel(data)
    if data{i+1} == -1
        fprintf(fid,'%s', data{i});
        break
    else
        fprintf(fid,'%s\n', data{i});
    end
end
fclose(fid);

end


function newstr = add_string (str,args,delimiter)

fun = @(x) strrep(x,delimiter,args);
newstr = cellfun(fun,str,'UniformOutput',false);

end