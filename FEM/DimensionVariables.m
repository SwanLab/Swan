classdef DimensionVariables < handle
    

    properties (Access = private)
        pdim
        mesh
    end
    
    methods (Static, Access = public)

        function obj = create(cParams)
            switch cParams.type
                case 'Scalar'
                    obj = DimensionScalar(cParams);
                case 'Vector'
                    obj = DimensionVector(cParams);
            end
        end
    end
    
%     methods (Access = public)
% 
%         function obj = DimensionVariables(cParams)
%             obj.init(cParams);
%         end
% 
%         function applyNdimfield(obj, num)
%             obj.ndimField = num;
%             obj.ndofPerElement = obj.nnodeElem*obj.ndimField;
%             obj.ndof           = obj.mesh.nnodes*obj.ndimField;
%         end
% 
%     end
%     
%     methods (Access = private)
% 
%         function obj = init(obj, cParams)
%             obj.mesh  = cParams.mesh;
%             obj.pdim  = cParams.pdim;
%         end        
%       
%     end
%     
end
