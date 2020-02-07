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
                
            'testPlotCylinderTetrahedra';
            'testPlotCylinderHexahedra';            
                
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
            

            };
        end
    end
end
