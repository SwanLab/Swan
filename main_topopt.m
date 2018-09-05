clc; clear variables; close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

% %% Test
% run('test_fem.m');
run('test_topopt.m');
% run('test_integration.m')
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
%         'GrippingTetrahedraCoarse_Case_1_1_1';
%         'GrippingTetrahedraCoarse_Case_2_1_1';
    %     'GrippingTetrahedraCoarse_Case_3_1_1';
    %     'GrippingTetrahedraCoarse_Case_4_1_1';
    %     'GrippingTetrahedraCoarse_Case_1_2_1';
    %     'GrippingTetrahedraCoarse_Case_2_2_1';
    %     'GrippingTetrahedraCoarse_Case_3_2_1';
    %     'GrippingTetrahedraCoarse_Case_4_2_1'
%         'CantileverQuadrilateral_Case_1_2_1';
%         'CantileverQuadrilateral_Case_1_2_2'
%         'CantileverQuadrilateral_Case_5_2_1'
%         'CantileverTriangle_Case_1_2_1'
%         'CantileverTriangle_Case_2_2_1'
%         'CantileverTriangle_Case_3_2_1'
%         'CantileverTriangle_Case_1_2_4'
%         'CantileverTriangle_Case_4_1_2'
%  'CantileverTriangle_Case_1_2_1'
%         'BridgeQuadrilateral_Case_5_1_1'
%         'BridgeQuadrilateral_Case_5_2_1'
%         'BridgeQuadrilateral_Allaire'
%         'BridgeQuadrilateral_Case_5_3_1'
%         'CantileverHexahedra_Case_1_1_1'
% 'CantileverHexahedra_Case_1_1_2'
%         'CantileverHexahedra_Case_5_1_1'
%         'CantileverHexahedra_Case_5_2_1'
    %     'CantileverHexahedra_Case_5_2_2'
    %     'CantileverHexahedra_Case_5_2_3'
%         'CantileverHexahedra_Case_5_1_2'
%         'CantileverHexahedra_Case_5_1_3'
%         'CantileverHexahedra_Case_5_1_4'
%         'CantileverHexahedra_Case_5_1_5'
    %     'CantileverHexahedra_Case_5_1_6'
    %     'CantileverHexahedra_Case_5_1_7'
%         'CantileverHexahedra_Case_5_1_8'
%         'CantileverHexahedra_Case_5_1_9'
%         'CantileverHexahedra_Case_5_1_10'
%         'SphereHexahedra_Test_2'
%         'SphereHexahedra_Test_4'
%         'SphereHexahedra_Test_8'
%         'SphereHexahedra_Test_16'
%         'SphereHexahedra_Test_32'
%     'SphereTetrahedra_Test_8'
% 'BridgeTetrahedraCoarse_Case_1_1_1'
% 'CantileverTetrahedra_Case_5_1_1'
% 'CantileverTetrahedraCoarse_Case_1_1_1'
% 'CantileverTetrahedra_Case_5_1_2'
% 'test_cantilever3'
% 'BridgeQuadrilateral_Case_5_2_1'
    };
for icases=1:size(filenames,1)
    clearvars -except filenames icases;
    close all;
    settings=Settings(filenames{icases});
    
    test = TopOpt_Problem(settings);
    test.preProcess;
    test.computeVariables;
    test.postProcess;
end