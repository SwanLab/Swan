classdef FE_FunctionFactory < handle
    
    methods (Access = public, Static)
        
        function FE_Function = create(mesh,dim,cParams)
            type = cParams.type;
            order = cParams.order;
            
            switch type
                case "Lagrange Simplicial"
                    FE_Function = FE_LagrangianFunction.create(mesh,dim,order);
                case "Raviart-Thomas"
                    FE_Function = FE_RaviartThomasFunction.create(mesh,dim,order);
                case "Nedelec"
                    FE_Function = FE_NedelecFunction.create(mesh,dim,order);
            end
        end
        
    end
end

