classdef GeometryType < handle
    
    methods (Static, Access = public)
        
        function g = compute(cParams)
            ndim  = size(cParams.coord,2);
            nGeom = ndim + cParams.kFace;
            switch nGeom
                case 1
                    g = 'Line';
                case 2
                    g = 'Surface';
                case 3
                    g = 'Volume';
            end
        end
        
    end
    
end