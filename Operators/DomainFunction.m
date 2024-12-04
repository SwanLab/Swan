classdef DomainFunction < BaseFunction

    properties (GetAccess = public, SetAccess = protected)
        operation
    end
    
    methods (Access = public)
        
        function obj = DomainFunction(cParams)
            obj.init(cParams)
        end
        
        function r = evaluate(obj,xV)
            r = obj.operation(xV);
        end

        function plot(obj)
            fD = obj.project('P1D',obj.mesh);
            fD.plot();
        end

        function fun = project(obj,target,mesh)
            s.mesh          = mesh;
            s.projectorType = target;
            proj = Projector.create(s);
            fun = proj.project(obj);
        end


    end

   methods (Access = private)

        function init(obj,cParams)
            obj.operation = cParams.operation;
            if isfield(cParams,'ndimf')
                obj.ndimf = cParams.ndimf;
            else
                obj.ndimf = 1;
            end

            if isfield(cParams,'mesh')
                obj.mesh = cParams.mesh;
            else
                obj.mesh = [];
            end
        end

    end
    
    


end