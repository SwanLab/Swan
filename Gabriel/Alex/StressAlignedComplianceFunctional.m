classdef StressAlignedComplianceFunctional < handle
    properties (Access = private)
        value0
        lastDesign
    end
    properties (Access = private)
        mesh
        filters
        compliance
        C
        dC
        stateProblem
        orientation
    end
    methods (Access = public)
        function obj = StressAlignedComplianceFunctional(cParams)
            obj.init(cParams);
        end
        function [J,dJ] = computeFunctionAndGradient(obj,x)
            obj.updateOrientation(x);
            xD     = x.obtainDomainFunction();
            xR     = obj.filterFields(xD);
            fields = [xR(:);{obj.orientation.getFunction()}];
            Cx     = obj.C(fields);
            dCx    = obj.dC(fields);
            dCx    = dCx(1:numel(xR));
            dxC    = ChainRule.compute(x,dCx);
            [J,dJ] = obj.compliance.computeFunctionAndGradient(Cx,dxC);
            dJ     = obj.filterFields(dJ);
            if isempty(obj.value0)
                obj.value0 = J;
            end
            J = J/obj.value0;
            for i = 1:numel(dJ)
                dJ{i}.setFValues(dJ{i}.fValues/obj.value0);
            end
        end
        function J0 = getInitialValue(obj)
            J0 = obj.value0;
        end
    end
    methods (Access = private)
        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.filters      = cParams.filters;
            obj.C            = cParams.C;
            obj.dC           = cParams.dC;
            obj.compliance   = cParams.complainceFromConstitutive;
            obj.stateProblem = cParams.stateProblem;
            obj.orientation  = cParams.orientation;
        end
        function updateOrientation(obj,x)
            design = x.fun.fValues;
            if isequal(design,obj.lastDesign)
                obj.orientation.update(obj.stateProblem.stressFun);
            end
            obj.lastDesign = design;
        end
        function xR = filterFields(obj,x)
            xR = cell(size(x));
            for i = 1:numel(x)
                xR{i} = obj.filters{i}.compute(x{i},2);
            end
        end
    end
    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Compliance';
        end
    end
end