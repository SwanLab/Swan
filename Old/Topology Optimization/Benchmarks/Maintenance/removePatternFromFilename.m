clear; close all; clc;

%% ***************** REMOVE PATTERN FROM FILENAME CASES **************** %%
% Remove a certain pattern from the name of a file.

folderpath = uigetdir;
list = updateList(folderpath);

pattern_to_remove = 'ORIOL_';

for i = 1:length(list)
    old_name = list(i).name;
    if contains(old_name,pattern_to_remove)
        new_name = getNewName(old_name,pattern_to_remove);
        new_name = fixIndex(new_name,folderpath);
        movefile(fullfile(folderpath,old_name),fullfile(folderpath,new_name));
        list2 = updateList(folderpath);
    end
end


function new_name = getNewName(old_name,pattern_to_remove)
    pos = strfind(old_name,pattern_to_remove);
    new_name = old_name(pos+length(pattern_to_remove):end);
end