function ComputingTopOpt

% Note: You can use FEMInputWriter to create benchmarking tests!
% Note: Use gid to create harder tests!
bool = 1;
if bool==1
    s.testName = 'test_micro3d';
    %s.testName = 'test_nullspace';
    t = TopOptComputer(s);
    t.compute();
    % With the following lines you obtain the result for the last iteration
    % (example: design variable with GiD. Test other results and also ParaView!)
    p1Params.fValues = t.computation.designVariable.value;
    p1Params.mesh    = t.computation.designVariable.mesh;
    Result           = P1Function(p1Params);
    unfitted_mesh = t.computation.designVariable.getUnfittedMesh;

    %for 2d and export
    %IMcond = unfitted_mesh.createInnerMeshGoodConditioning();
    %height = 0.16;
    %EM = IMcond.provideExtrudedMesh(height);
    %EM.exportSTL();

    %for 3d - doesnt work
    %unfitted_mesh.exportSTL();

    %c.type = 'GiD';
    %c.filename = [s.testName,'_LastIter'];
    %Result.print(c);
end

dir_init = strcat(cd,'/','Output/test_micro3d');
dir_server = 'P:/03_Foerderprojekte/Bayern/BStMWi-CC40-MAI_rapidSkelett/03_Simulation/2_Runs/SWAN/autocopyCalib/';
folder_server = 'OScomp_AF_aJ_2.5_aG_015_MINRES_SphereRad_04_v045_mesh005_P1_msmooth_it8';
path_server = append(dir_server, folder_server);
if exist(path_server,'dir') == 0, mkdir(path_server); end
myFiles = dir(fullfile(dir_init,'*.vtu'));
for k = 1:(length(myFiles))
  baseFileName = myFiles(k).name;
  fullFileName = fullfile(dir_init, baseFileName);
  copyfile(fullFileName,path_server);
  delete(fullFileName);
end
