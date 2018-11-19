%% INTEGRATION TEST =======================================================

clear; close all;

fprintf('Running TopOpt tests...\n')

%% Test Declaration -------------------------------------------------------

tests_integration = {'test_circle_triangle','test_circle_quadrilateral','test_sphere_tetrahedra','test_sphere_hexahedra','test_cylinder_tetrahedra','test_cylinder_hexahedra'};
% tests_integration = {'test_cylinder_tetrahedra','test_cylinder_hexahedra'};

%% Run Integration Opt Tests ----------------------------------------------
for i = 1:length(tests_integration)
    clearvars -except tests_integration i
    %     tic
    file_name = tests_integration{i};
    file_name_in = strcat('./Input/',tests_integration{i});
    settings = Settings(file_name_in);
    load_file = strcat('./tests/',file_name);
    load(load_file)
    
    obj = TopOpt_Problem(settings);
    obj.preProcess;
    x = obj.x;
    
    mesh_boundary = Mesh_Unfitted.create('BOUNDARY',obj.mesh,Interpolation.create(obj.mesh,'LINEAR'));
    mesh_boundary.computeMesh(x);
    
    switch obj.mesh.pdim
        case '2D'
            A = mesh_boundary.computePerimeter;
        case '3D'
            A = mesh_boundary.computeSurface;
    end
    
    errorSurf = A/A0 - 1;
    if abs(errorSurf) < 4e-2
        cprintf('green',strcat(file_name,' PASSED.  Surface Error: ',num2str(errorSurf),'\n'));
    else
        cprintf('err',strcat(file_name,' FAILED. Surface Error: ',num2str(errorSurf),'\n'));
    end
    
    mesh_interior = Mesh_Unfitted.create('INTERIOR',obj.mesh,Interpolation.create(obj.mesh,'LINEAR'));
    mesh_interior.computeMesh(x);
    
    switch obj.mesh.pdim
        case '2D'
            V = mesh_interior.computeSurface;
        case '3D'
            V = mesh_interior.computeVolume;
    end
    
    errorVol= V/V0 - 1;
    if abs(errorVol) < 6e-2
        cprintf('green',strcat(file_name,' PASSED.  Volume Error: ',num2str(errorVol),'\n'));
    else
        cprintf('err',strcat(file_name,' FAILED. Volume Error: ',num2str(errorVol),'\n'));
    end
    
    %     toc
    clear settings
end

fprintf('\nTopOpt tests completed.\n')
fprintf('\n-------------------------------------------\n\n')