classdef FE_FunctionFactory < handle
    
    methods (Access = public, Static)
        
        function FE_Function = create(mesh,dim,cParams)
            type = cParams.type;
            order = cParams.order;
            finiteElement = cParams.finiteElement;
            
            switch type
                case "Lagrange Simplicial"
                    FE_Function = FE_LagrangianFunction.create(mesh,dim,order,finiteElement);
                case "Raviart-Thomas"
                    FE_Function = FE_RaviartThomasFunction.create(mesh,dim,order,finiteElement);
                case "Nedelec"
                    FE_Function = FE_NedelecFunction.create(mesh,dim,order,finiteElement);
            end
        end
        
        function FE_Function = createWithFValues(s)
            type = s.type;
            
            switch type
                case "Lagrange Simplicial"
                    FE_Function = FE_LagrangianFunction(s);
                case "Raviart-Thomas"
                    FE_Function = FE_RaviartThomasFunction(s);
                case "Nedelec"
                    FE_Function = FE_NedelecFunction(s);
            end
        end
        
    end
end

