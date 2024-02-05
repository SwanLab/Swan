classdef LagrangeElement < handle
    
    methods (Access = public, Static)
        
        function element = create(elementType,k,d)
            switch elementType
                case "SIMPLICIAL"
                    if d==1
                        element = LagrangeSimplicial1D(k);
                    elseif d==2
                        element = LagrangeSimplicial2D(k);
                    elseif d==3
                        element = LagrangeSimplicial3D(k);
                    end
                case "TENSOR PRODUCT"
                    if d==1
                        element = LagrangeTensorProduct1D(k);   
                    elseif d==2
                        element = LagrangeTensorProduct2D(k);
                    elseif d==3
                        element = LagrangeTensorProduct3D(k);
                    end
                case "PRISMATIC"
                    if d==1
                        element = LagrangePrismatic1D(k);   
                    elseif d==2
                        element = LagrangePrismatic2D(k);
                    elseif d==3
                        element = LagrangePrismatic3D(k);
                    end
            end
        end
        
    end
    
    methods (Access = public, Abstract)
        
    end
    
end