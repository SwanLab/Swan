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
          
                                   
            
            'testSmoothRectangleTriangle'                          
            'testSmoothRectangleQuadrilateral'
                
                            
                        
            
            'testPlotCircleQuadrilateral'   
            'testPlotCircleTriangle'                                
            
          
            
            'testRectangleTriangle'
            'testRectangleQuadrilateral'           
            
            'testCircumferenceTriangle'
            'testCircumferenceQuadrilateral'
            
 
            'testPlotCylinderTetrahedra';
            'testPlotCylinderHexahedra';      
                       
            
            
            'testPlotSphereTetrahedra';
            'testPlotSphereHexahedra';                                
        
            
             
    

            };
        end
    end
end
