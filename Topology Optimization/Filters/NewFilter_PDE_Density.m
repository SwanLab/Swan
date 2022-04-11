classdef NewFilter_PDE_Density < handle
    
    properties (Access = private)
        epsilon
        diffReacProb
    end

    methods (Access = public)
        
        function obj = NewFilter_PDE_Density(cParams)
            obj.init(cParams);
            obj.epsilon = cParams.mesh.computeMeanCellSize();
        end
        
        function preProcess(obj)
            obj.diffReacProb.setEpsilon(obj.epsilon);
        end

        function RHS = integrate_L2_function_with_shape_function(obj,x)
            RHS = obj.diffReacProb.getM()*x;
        end

        function obj = updateEpsilon(obj,epsilon)
            if obj.hasEpsilonChanged(epsilon)
                obj.epsilon = epsilon;
                obj.diffReacProb.setEpsilon(epsilon);
                obj.computeLHS();
            end
        end

        function x_reg = getP1fromP1(obj,x)
            RHS = obj.integrate_L2_function_with_shape_function(x);
            x_reg = obj.solve_filter(RHS);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.createDiffReacProblem(cParams);
        end

        function createDiffReacProblem(obj,cParams)
            s = cParams.femSettings;
            s.mesh = cParams.mesh;
            switch s.scale
                case 'MACRO'
                    obj.diffReacProb = DiffReactProblem(s);
                case 'MICRO'
                    obj.diffReacProb = DiffReactProblemMicro(s);
            end
        end

        function x_reg = solve_filter(obj,RHS)
            obj.diffReacProb.computeVariables(RHS);
            x_reg = obj.diffReacProb.variables.x;
        end

        function itHas = hasEpsilonChanged(obj,eps)
            if isempty(obj.epsilon)
                obj.epsilon = 0;
            end
            var = abs(eps - obj.epsilon)/eps;
            itHas = var > 1e-15;
        end

    end

end