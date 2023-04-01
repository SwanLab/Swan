classdef DesignVolumen < handle 
    properties (Access = public)
        volumenFraction
        volumen
        derivedVolumen
    end
    properties (Access = private)
        elementNumberX
        elementNumberY
        filterParameters
        projectedField
        derivedProjectedField       
    end
    methods (Access = public)
        function obj = DesignVolumen(cParams)
            obj.inputData(cParams);
        end
        function computeVolumen(obj)
            obj.volumen  = sum(sum(obj.projectedField))/(obj.volumenFraction*obj.elementNumberX*obj.elementNumberY)-1;
        end 
        function deriveVolumen(obj)
            obj.derivedVolumen  = 1/(obj.volumenFraction*obj.elementNumberX*obj.elementNumberY)*ones(obj.elementNumberY,obj.elementNumberX);
            obj.filterDerivedVolumen();
        end 
        function computeVolumenFraction(obj,D,I)
           obj.volumenFraction = obj.filterParameters.volumenFraction*sum(D.designField.projectedField(:))/sum(I.designField.projectedField(:));
        end 

    end
    methods (Access = private)
        function inputData(obj,cParams)
            obj.elementNumberX = cParams.mesh.elementNumberX;
            obj.elementNumberY =cParams.mesh.elementNumberY;
            obj.filterParameters =cParams.filterParameters;
            obj.projectedField =cParams.projectedField;
            obj.derivedProjectedField = cParams.derivedProjectedField;
            obj.volumenFraction = cParams.volumenFraction;
        end

        function filterDerivedVolumen(obj)
            obj.derivedVolumen(:)  =  obj.filterParameters.H*(obj.derivedVolumen(:).*obj.derivedProjectedField(:)./obj.filterParameters.Hs);
        end 
    end
end