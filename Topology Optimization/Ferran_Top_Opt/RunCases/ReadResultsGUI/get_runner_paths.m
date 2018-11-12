function user_paths = get_runner_paths

Runner = char(java.net.InetAddress.getLocalHost.getHostName);
switch Runner
    case 'FERRY58-PC'
        user_paths.default_folder = 'E:\TFM';
        user_paths.temp_folder = 'C:\Users\ferry\Downloads';
        user_paths.default_stl_path = 'E:\TFM';
        user_paths.macro_source = 'E:\Dropbox\Ferran\RunCases\ReadResultsGUI\Macros.tcl';
        user_paths.macro_dest = 'C:\Users\ferry\AppData\Roaming\GiD\Macros.tcl';
        user_paths.gid_path = 'gid';
        user_paths.dxf_path = 'C:/Program Files/GiD/GiD 13.0.1/templates';
        working_path = 'E:\Dropbox\Ferran\CodeTopOpt';
        copyfile(user_paths.macro_dest,user_paths.macro_source); % to keep the code updated in dropbox

    case 'alex'
        working_path = '/home/alex/Dropbox/Ferran/CodeTopOpt';
        user_paths.default_folder = '/home/alex/Desktop/Results';
        user_paths.temp_folder = '/home/alex/Downloads';
        user_paths.default_stl_path = '/home/alex/Dropbox/Ferran/STL models';
        user_paths.macro_dest = '/home/alex/Dropbox/Ferran/RunCases/ReadResultsGUI/Macros.tcl';
        user_paths.gid_path = '/opt/GiDx64/13.1.2d/gid';
        user_paths.dxf_path = '/opt/GiDx64/13.1.2d/templates';


    case 'GOKU'
        working_path = '/home/aferrer/Dropbox/Ferran/CodeTopOpt';

    otherwise
        edit get_runner_paths.m
        error('Runner not detected. Please, add runner: %s',Runner);
end

funpath = mfilename('fullpath');
funname = mfilename;
user_paths.mfilepath = funpath(1:end-length(funname));
addpath(genpath(working_path));

% EXAMPLE
% user_paths.default_folder = []; % default results path
% user_paths.temp_folder = []; % used for copying temp .png, .res, .msh files
% user_paths.default_stl_path = []; % default path for saving STL models
% user_paths.macro_dest = []; % folder of the TCL code
% user_paths.gid_path = []; % command to call gid
% user_paths.dxf_path = []; % exporting STL
end