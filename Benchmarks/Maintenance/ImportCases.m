clear; close all; clc;

%% *************************** IMPORT CASES **************************** %%
% Move cases from an origin folder to a destination folder considering case
% indexes.

origin_folderpath = uigetdir;
destination_superfolderpath = uigetdir;
list = updateList(origin_folderpath);

for i = 1:length(list)
    old_name = list(i).name;
    destination_folderpath = findDestinationFolder(old_name,destination_superfolderpath);
    new_name = fixIndex(old_name,destination_folderpath);
    if ~exist(destination_folderpath,'dir')
        mkdir(destination_folderpath);
    end
    movefile(fullfile(origin_folderpath,old_name),fullfile(destination_folderpath,new_name));
end


function destination_folderpath = findDestinationFolder(name,destination_superfolderpath)
    
    if contains(name,'Bridge','IgnoreCase',true)
        destination_folderpath = fullfile(destination_superfolderpath,'Bridge');
    elseif contains(name,'Cantilever','IgnoreCase',true)
        destination_folderpath = fullfile(destination_superfolderpath,'Cantilever');
    elseif contains(name,'Throne','IgnoreCase',true)
        destination_folderpath = fullfile(destination_superfolderpath,'Throne');
    elseif contains(name,'Chair','IgnoreCase',true)
        destination_folderpath = fullfile(destination_superfolderpath,'Chair');
    else
        error('NOT FOUND')
    end
end