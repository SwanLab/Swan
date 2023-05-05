classdef FiniteElementFactory < handle
    
    methods (Access = public, Static)
        
        function finiteElement = create(cParams)
            type = cParams.type;
            order = cParams.order;
            dim = cParams.dim;
            
            switch type
                case "Lagrange Simplicial"
                    finiteElement = LagrangeElement.create("SIMPLICIAL",order,dim);
                case "Raviart-Thomas"
                    finiteElement = RaviartThomasElement.create(dim);
                case "Nedelec"
                    finiteElement = NedelecElement.create(dim);
            end
        end
        
    end
end

