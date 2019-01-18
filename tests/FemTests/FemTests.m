classdef FemTests < testRunner
    
    
    properties (Access = protected)
        FieldOfStudy = 'FEM'
        tests
    end
    
    methods (Access = public)
        function obj = FemTests()
            obj@testRunner();
        end
    end
    
    methods (Access = protected)
        function loadTests(obj)
            obj.tests = {...  
                'test2dQuad';                
                'testStiffnessMatrixGenerator';
                'test2dTriangle';
                'test2dStokes_triangle';
                'test2dMicro';
                'test3dTetrahedra';
                'test3dHexahedra'
                };

        end
    end
    
end

