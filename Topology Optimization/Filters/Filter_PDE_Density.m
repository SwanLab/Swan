classdef Filter_PDE_Density < Filter_PDE
    
    methods (Access = public)
        
        function obj = Filter_PDE_Density(cParams)
            obj.init(cParams);
            obj.epsilon = cParams.mesh.computeMeanCellSize();
            %obj.epsilon = 0;
        end
        
        function preProcess(obj)
            preProcess@Filter(obj)
            obj.Anodal2Gauss = obj.computeA();
%             obj.diffReacProb.setEpsilon(obj.epsilon);
            obj.computeLHS(obj.epsilon);
        end
        
        function RHS = integrate_L2_function_with_shape_function(obj,x)
            dim = obj.computeDimensions();
            M = obj.computeMassMatrix(dim);
            RHS = M*x;
        end
        
    end

    methods (Access = private)

        function dim = computeDimensions(obj)
            s.ngaus = [];
            s.mesh  = obj.mesh;
            s.pdim  = '1D';
            d       = DimensionVariables(s);
            d.compute();
            dim = d;
        end
        
        function M = computeMassMatrix(obj, dim)
            s.type         = 'MassMatrix';
            s.quadType     = 'QUADRATICMASS';
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.dim          = dim;
            LHS = LHSintegrator.create(s);
            M = LHS.compute();
        end

    end
    
end