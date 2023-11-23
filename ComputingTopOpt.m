function ComputingTopOpt

% Note: You can use FEMInputWriter to create benchmarking tests!
% Note: Use gid to create harder tests!
bool = 1;
if bool==1
    s.testName = 'test_micro3d';
    t = TopOptComputer(s);
    t.compute();
    % With the following lines you obtain the result for the last iteration
    % (example: design variable with GiD. Test other results and also ParaView!)
    p1Params.fValues = t.computation.designVariable.value;
    p1Params.mesh    = t.computation.designVariable.mesh;
    Result           = P1Function(p1Params);
    c.type = 'GiD';
    c.filename = [s.testName,'_LastIter'];
    Result.print(c);
end

dir_init = strcat(cd,'/','Output/test_micro3d');
dir_server = 'P:/03_Foerderprojekte/Bayern/BStMWi-CC40-MAI_rapidSkelett/03_Simulation/2_Runs/SWAN/autocopy/';
folder_server = 'auto_[1 0.125 0.125 0.3 0.125 0.3]_aJ_2.5_aG_015_MINRES_SphereRad_06_v027_mesh0025';
path_server = append(dir_server, folder_server);
if exist(path_server,'dir') == 0, mkdir(path_server); end
myFiles = dir(fullfile(dir_init,'*.vtu'));
for k = 1:(length(myFiles))
  baseFileName = myFiles(k).name;
  fullFileName = fullfile(dir_init, baseFileName);
  copyfile(fullFileName,path_server);
  delete(fullFileName);
end
