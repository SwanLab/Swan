classdef UnfittedIntegrationTests < testRunner
    
    properties (Access = protected)
        FieldOfStudy = 'Unfitted integration'
        tests
    end
    
    methods (Access = public)
        
        function  obj = UnfittedIntegrationTests()
            obj@testRunner();
        end
        
    end
    
    methods (Access = protected)
        
        function loadTests(obj)
            
            
            obj.tests = {...

               

                      
             'testSurfaceSphereTetrahedra';                       

                     
            'testPerimeterRectangleTriangle';                
            'testPerimeterRectangleQuadrilateral';                      
            
            'testPerimeterCircleTriangle';
            'testPerimeterCircleQuadrilateral';
            
            'testAreaCircleTriangle';            
            'testAreaCircleQuadrilateral';               
                
           'testSurfaceSphereHexahedra';   
             
           'testSurfaceCylinderTetrahedra';
           'testSurfaceCylinderHexahedra';                   
                
           'testVolumeSphereTetrahedra';
           'testVolumeSphereHexahedra';  
           
           'testVolumeCylinderTetrahedra';
           'testVolumeCylinderHexahedra';                  

            
                          

    
    
            };
        end
        
    end
    
end

