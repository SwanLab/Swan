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
                
            
            'testPerimeterRectangleQuadrilateral';                      
            
             'testPerimeterCircleQuadrilateral';
            
            'testAreaCircleQuadrilateral';       
            
           'testSurfaceSphereHexahedra';   
           'testSurfaceSphereTetrahedra';                       

           'testVolumeSphereTetrahedra';
           'testVolumeSphereHexahedra';  
                   
            

            'testPerimeterRectangleTriangle';                            
            'testPerimeterCircleTriangle';            
       
            'testAreaCircleTriangle';                 

            

           'testSurfaceCylinderTetrahedra';            
           'testSurfaceCylinderHexahedra';   
           
           'testVolumeCylinderTetrahedra';
           'testVolumeCylinderHexahedra';              
                                  

     
    

                                  
                          
       
    
            };
        end
        
    end
    
end

