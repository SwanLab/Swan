classdef PlottingTests < testRunner
    properties (Access = protected)
        FieldOfStudy = 'Plotting tests'
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
                %                 'testRectangleTriangle'
                %                 'testRectangleQuadrilateral'
                
                'testPlotCylinderTetrahedra';
                'testPlotCylinderHexahedra';
                
                'testPlotCircleTriangle'
                'testPlotCircleQuadrilateral'
                
                'testPlotSphereTetrahedra';
                'testPlotSphereHexahedra';
                };
        end
    end
end
