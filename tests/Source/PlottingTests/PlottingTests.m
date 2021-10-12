classdef PlottingTests < testRunner
    properties (Access = protected)
        FieldOfStudy = 'Plotting'
        tests
    end
    
    methods (Access = public)
        function  obj = PlottingTests()
            obj@testRunner();
            close all;
        end
    end
    
    methods (Access = protected)
        function loadTests(obj)
            obj.tests = {... 
            %% coords, levelset especificats (PlottingToyUnfittedExample)
%             'testTriangleToyUntittedExample'
%             'testQuadToyUntittedExample'
            %% meshContour, boundary/int... (testPlotting)
%             'testCircumferenceQuadrilateral'
%             'testCircumferenceTriangle'
%             'testPlotCircleTriangle'
%             'testPlotCircleQuadrilateral'
% 
%             'testPlotSphereTetrahedra';
%              'testPlotSphereHexahedra';
            
            %% meshContour, boundary/int... (testPlotting_Composite)
%             'testRectangleTriangle'
%             'testRectangleQuadrilateral'
%             'testSmoothRectangleTriangle'
%             'testSmoothRectangleQuadrilateral'
% 
%              'testPlotCylinderTetrahedra';
%              'testPlotCylinderHexahedra';
             %% bolet
             'testPlotLargeCylinderTethaedra';
            };
        end
    end
end
