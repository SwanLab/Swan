classdef OneTetrahedronRandCoordOrderedConnecCase7 <  VectorizedTriangulationTest
    
    properties (Access = private)
       testName = 'OneTetrahedronRandCoordOrderedConnecCase7';        
    end
    
    methods (Access = protected)
        
        function init(obj)
            obj.coord = rand(4,3);
%             obj.coord =  [0.5938    0.2836    0.9047;
%                          0.2827    0.5508    0.1310;
%                          0.1552    0.8709    0.8337;
%                          0.0007    0.0423    0.8005;];

         %   obj.coord =  [0 0 0;0 1 0;1 0 0;0 0 1];
            obj.connec = [1 2 3 4];
            obj.boundaryConnec = [1 2 3;1 2 4;1 3 4;2 3 4];
            obj.levelSet = [-7.8496;-9.7731;-8.3404;8.3622];
            obj.connecBcutMesh = [5 6 7];                        
        end
        
        
        
    end
    
end