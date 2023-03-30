classdef DesignField < handle
    properties (Access = public)
        field
        projectedField
        filteredField
        derivedProjectedField
        derivedFilteredField
    end
    properties (Access = private)
        filterParameters
        projectorParameters
        mesh
    end

    methods (Access = public)
        function obj = DesignField(cParams)
            obj.inputData(cParams)
            obj.computeWeigth()
        end
        function filter(obj)
            s.filterParameters = obj.filterParameters;
            s.field = obj.field;
            B = FilterComputer(s);
            B.compute();
            obj.filteredField = B.filteredField;
        end
        function project(obj)
            if isempty(obj.filteredField) == 1 
                warning("Field should be filtered before projecting")
            end
            s.beta = obj.projectorParameters.beta;
            s.eta = obj.projectorParameters.eta;
            s.filteredField =obj.filteredField;
            B = FieldProjector(s);
            B.compute();
            obj.projectedField = B.projectedField;
        end
        function deriveProjectedField(obj)
            s.beta = obj.projectorParameters.beta;
            s.filteredField = obj.filteredField;
            s.eta =obj.projectorParameters.eta;
            B = ProjectedFieldFilteredFieldDerivator(s);
            B.compute();
            obj.derivedProjectedField = B.derivatedProjectedField;
        end
        function deriveFilteredField(obj)
            s.H = obj.filterParameters.H;
            s.Hs = obj.filterParameters.Hs;
            B = FilteredFieldFieldDerivator(s);
            B.compute();
            obj.derivedFilteredField = B.derivedFilteredField;
        end
    end
    methods (Access = private)
        function inputData(obj,cParams)
            obj.field = cParams.field;
            obj.projectorParameters =  cParams.projectorParameters;
            obj.filterParameters =  cParams.filterParameters;
            obj.mesh = cParams.mesh;
        end
        function computeWeigth(obj)
            s.elementNumberX =obj.mesh.elementNumberX;
            s.elementNumberY = obj.mesh.elementNumberY;
            s.minimunInfluenceRadios = obj.filterParameters.minimunInfluenceRadios;
            B = weightFilterComputer(s);
            B.compute();
            obj.filterParameters.H = B.H;
            obj.filterParameters.Hs = B.Hs;
        end  
    end
end