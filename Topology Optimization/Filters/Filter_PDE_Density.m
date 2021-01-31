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
            obj.diffReacProb.setEpsilon(obj.epsilon);
            obj.computeLHS();
        end        
        
        function RHS = integrate_L2_function_with_shape_function(obj,x)
            RHS = obj.diffReacProb.element.M*x;
        end
        
    end
    
end