classdef Tetrahedra_Linear<Interpolation
    
    methods (Access = public)
        
        function obj = Tetrahedra_Linear(cParams)
            obj.init(cParams);
            obj.computeParams()
            obj.computeCases();
            obj.computeCoordAndConnec();
        end
        
        function computeShapeDeriv(obj,posgp)
            obj.shape = [];
            s = posgp(1,:);
            t = posgp(2,:);
            u = posgp(3,:);
            obj.shape =[(ones(1,size(posgp,2))-t-s-u);
                s;
                t;
                u];
            obj.deriv = repmat([-1 1 0 0
                -1 0 1 0
                -1 0 0 1],1,1,size(posgp,2));            
        end
        
    end
    
    methods (Access = private)
        
        function computeParams(obj)
            obj.type  = obj.mesh.geometryType;
            obj.ndime = 3;
            obj.nnode = 4;
            obj.pos_nodes = [0 0 0;
                1 0 0;
                0 1 0;
                0 0 1];
            obj.isoDv = 1/6;
        end
        
        function computeCases(obj)
            obj.iteration = [1 1 1 2 2 3;
                2 3 4 3 4 4];
            obj.cases(:,:,1) = [7 6 5 1
                5 6 7 4
                5 6 4 2
                4 3 2 6;
                zeros(2,4)];
            obj.cases(:,:,2) = [5 6 7 2
                7 6 5 3
                5 7 3 1
                7 3 1 4;
                zeros(2,4)];
            obj.cases(:,:,3) = [7 6 5 3
                5 6 7 2
                2 7 5 1
                2 7 1 4;
                zeros(2,4)];
            obj.cases(:,:,4) = [5 6 7 4
                7 6 5 3
                6 5 3 1
                6 1 3 2;
                zeros(2,4)];
            obj.cases(:,:,5) = [2 5 7 3
                5 8 7 3
                5 6 8 3
                7 8 1 4
                7 6 5 1
                7 8 6 1];
            obj.cases(:,:,6) = [6 8 3 4
                6 8 7 3
                7 5 6 3
                6 5 2 1
                8 5 7 2
                8 6 5 2];
            obj.cases(:,:,7) = [2 8 6 4
                5 8 6 2
                5 7 8 2
                6 3 5 1
                8 7 5 3
                6 8 5 3];
            obj.main_loop = [4 4];
            obj.extra_cases = [5,6,7];
            obj.selectcases = [1 0 0;
                2 0 0;
                3 6 0;
                4 7 0;
                0 5 0;
                0 7 4;
                0 6 3
                0 0 2
                0 0 1];
        end
    end
    
end
