%% To-do
% - GiDImageCapturer -> some variables should be "user variables", defined
%   once in a centralized file

% swanPath = '/home/ton/Github/Swan/';
% gidPath  = '/home/ton/GiDx64/gid-16.1.2d/';
% pathTcl  = [swanPath,'PostProcess/STL/'];
% 
% tclFile = [pathTcl,'CreateSurfaceMeshFile.tcl" '];
% command = [gidPath,'gid_offscreen -offscreen -t "source ',tclFile];
% templateText = fileread('CreateSurfaceMeshFile_Template.tcl');
% 
% % set input_post_res "/home/ton/test_micro23.flavia.res"
% % set output_gid_project_name "/home/ton/test_micro_project.gid"
% % set mesh_element_size "0.0707107"
% % set mesh_name "hmmmm"
% 
% fid = fopen('CreateSurfaceMeshFile.tcl', 'w');
% fprintf(fid, 'New message in new line\n');
% fclose(fid);
% 
% system(command);


s.filename        = 'hellothere';
s.gidProjectPath  = '/home/ton/test_micro_project.gid';
s.meshElementSize = '0.0707107';
s.meshFileName    = 'hmmmm22';
swanPath = '/home/ton/Github/Swan/';
gidPath  = '/home/ton/GiDx64/gid-16.1.2d/';