%% To-do
% - GiDImageCapturer -> some variables should be "user variables", defined
%   once in a centralized file

swanPath = '/home/ton/Github/Swan/';
gidPath  = '/home/ton/GiDx64/gid-16.1.2d/';
pathTcl  = [swanPath,'PostProcess/STL/'];

tclFile = [pathTcl,'CreateSurfaceMeshFile.tcl" '];
command = [gidPath,'gid_offscreen -offscreen -t "source ',tclFile];
system(command);