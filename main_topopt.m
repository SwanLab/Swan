clc; clear variables; close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

%% Test
run('test_fem.m');
run('test_topopt.m');
clear variables;
%% Main
filenames={%'GrippingTriangleCoarse_Case_1_1_1';
    %     'GrippingTriangleCoarse_Case_2_1_1';
    %     'GrippingTriangleCoarse_Case_3_1_1';
    %     'GrippingTriangleCoarse_Case_4_1_1';
    %     'GrippingTriangleFine_Case_1_1_1';
    %     'GrippingTriangleFine_Case_2_1_1';
    %     'GrippingTriangleFine_Case_3_1_1';
    %     'GrippingTriangleFine_Case_4_1_1';
    %     'GrippingQuadCoarse_Case_1_1_1';
    %     'GrippingQuadCoarse_Case_2_1_1';
    %     'GrippingQuadCoarse_Case_1_1_1';
    %     'GrippingQuadCoarse_Case_4_1_1';
    %     'GrippingQuadFine_Case_1_1_1';
    %     'GrippingQuadFine_Case_2_1_1';
    %     'GrippingQuadFine_Case_3_1_1';
    %     'GrippingQuadFine_Case_4_1_1';
    %     'GrippingTetrahedraCoarse_Case_1_1_1';
    %     'GrippingTetrahedraCoarse_Case_2_1_1';
    %     'GrippingTetrahedraCoarse_Case_3_1_1';
    %     'GrippingTetrahedraCoarse_Case_4_1_1';
    %     'GrippingTetrahedraCoarse_Case_1_2_1';
    %     'GrippingTetrahedraCoarse_Case_2_2_1';
    %     'GrippingTetrahedraCoarse_Case_3_2_1';
    %     'GrippingTetrahedraCoarse_Case_4_2_1'
    %     'CantileverQuadrilateral_Case_1_2_1';
    %     'CantileverQuadrilateral_Case_1_2_2'
    %     'CantileverQuadrilateral_Case_5_2_1'
    %     'CantileverTriangle_Case_1_2_1'
    %     'CantileverTriangle_Case_1_2_4'
    %     'BridgeQuadrilateral_Case_5_1_1'
    %     'BridgeQuadrilateral_Case_5_1_2'
    %     'BridgeQuadrilateral_Case_5_1_3'
    %     'CantileverTriangle_Case_4_1_2'
%         'BridgeQuadrilateral_Case_5_2_1'
    %     'BridgeQuadrilateral_Case_5_3_1'
    %     'CantileverHexahedra_Case_1_1_1'
    %     'CantileverHexahedra_Case_5_1_1'
    %     'CantileverHexahedra_Case_5_2_1'
    %     'CantileverHexahedra_Case_5_2_2'
    %     'CantileverHexahedra_Case_5_2_3'
    %     'CantileverHexahedra_Case_5_1_2'
%         'CantileverHexahedra_Case_5_1_3'
    %     'CantileverHexahedra_Case_5_1_4'
    %     'CantileverHexahedra_Case_5_1_5'
    %     'CantileverHexahedra_Case_5_1_6'
    %     'CantileverHexahedra_Case_5_1_7'
    %     'CantileverHexahedra_Case_5_1_8'
    %     'CantileverHexahedra_Case_5_1_9'
    %     'CantileverHexahedra_Case_5_1_10'
    %     'SphereHexahedra_Test_Case_5_2'
    %     'SphereHexahedra_Test_Case_5_4'
    %     'SphereHexahedra_Test_Case_5_8'
    %     'SphereHexahedra_Test_Case_5_16'
    %     'SphereHexahedra_Test_Case_5_32'
%     'test_bridge2'
    };
for icases=1:size(filenames,1)
    clearvars -except filenames icases;
    close all;
    settings=Settings(filenames{icases});
    % --------------------------- !! DELETE !! ----------------------------
    if contains(lower(filenames{icases}),'case_5') || contains(lower(filenames{icases}),'test_bridge2')
        if contains(lower(filenames{icases}),'test_bridge2')
            dim = [60 20 20]; div = [24 8 8];
            [A1,b1,A0,b0] = conversionTensors3D(settings.filename,dim,div);
        elseif ~contains(lower(filenames{icases}),'hexa')
            dim = [2 1]; div = [120 60];
            [A1,b1,A0,b0] = conversionTensors(settings.filename,dim,div);
        else
            if contains(lower(filenames{icases}),'cantilever')
                dim = [60 20 20]; div = [48 16 16];
                % dim = [60 20 20]; div = [96 32 32];
                [A1,b1,A0,b0] = conversionTensors3D(settings.filename,dim,div);
            elseif contains(lower(filenames{icases}),'sphere')
                dim = [2 2 2];
                if contains(filenames{icases},'32')
                    div = [32 32 32];
                elseif contains(filenames{icases},'16')
                    div = [16 16 16];
                elseif contains(filenames{icases},'8')
                    div = [8 8 8];
                elseif contains(filenames{icases},'4')
                    div = [4 4 4];
                elseif contains(filenames{icases},'2')
                    div = [2 2 2];
                else
                    dim = [1 1 1]; div = [25 25 25];
                end
                [A1,b1,A0,b0] = conversionTensors3D(settings.filename,dim,div);
            end
        end
        save(fullfile(pwd,'Allaire_ShapeOpt','conversion'),'A0','A1','b0','b1','dim','div');
    end
    % ---------------------------------------------------------------------
    test = TopOpt_Problem(settings);
    test.preProcess;
    test.computeVariables;
    test.postProcess;
end