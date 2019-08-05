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
                
            'testPlotCircleQuadrilateral'
            
            'testSmoothRectangleTriangle'
            'testSmoothRectangleQuadrilateral'
            
            'testRectangleTriangle'
            'testRectangleQuadrilateral'
            
            'testCircumferenceTriangle'
            'testCircumferenceQuadrilateral'
            
            'testPlotCircleTriangle'
            
            'testPlotSphereTetrahedra';
            'testPlotSphereHexahedra';
            
            'testPlotCylinderTetrahedra';
            'testPlotCylinderHexahedra';
            };
        end
    end
end
