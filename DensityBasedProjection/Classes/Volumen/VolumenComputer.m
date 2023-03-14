classdef VolumenComputer < handle 
    properties (Access = public)
        volumen
    end
    properties (Access = private)
        elementNumberX
        elementNumberY
        filterParameters
        projectedField
        derivedProjectedField       
    end
    methods (Access = public)
        function obj = VolumenComputer(cParams)
            obj.inputData(cParams);
        end

        function compute(obj)
            obj.computeVolumen();
            obj.deriveVolumen();
            obj.filterDerivedVolumen();
        end
    end
    methods (Access = private)
        function inputData(obj,cParams)
            obj.elementNumberX = cParams.mesh.elementNumberX;
            obj.elementNumberY =cParams.mesh.elementNumberY;
            obj.filterParameters =cParams.filterParameters;
            obj.projectedField =cParams.projectedField;
            obj.derivedProjectedField = cParams.derivedProjectedField;
            obj.volumen = cParams.volumen;
        end
        function computeVolumen(obj)
            obj.volumen.value        = sum(sum(obj.projectedField.D))/(obj.volumen.volfracD*obj.elementNumberX*obj.elementNumberY)-1;
        end 
        function deriveVolumen(obj)
            obj.volumen.derivated   = 1/(obj.volumen.volfracD*obj.elementNumberX*obj.elementNumberY)*ones(obj.elementNumberY,obj.elementNumberX);
        end
        function filterDerivedVolumen(obj)
            obj.volumen.derivatedFiltered(:)  = obj.filterParameters.H*(obj.volumen.derivated(:).*obj.derivedProjectedField.D(:)./obj.filterParameters.Hs);
        end 
    end
end