classdef Integrator_Composite < Integrator
    
    properties (GetAccess = public, SetAccess = private)
        integrators
        nInt
    end
    
    properties (Access = private)
       RHScells
       RHSsubcells
       boxFaceToGlobal
    end
    
    methods (Access = public)
        
        function obj = Integrator_Composite(cParams)
            obj.init(cParams);
            if isfield(cParams,'boxFaceToGlobal')
                obj.boxFaceToGlobal = cParams.boxFaceToGlobal;
            end
            obj.createIntegrators(cParams);
        end
        
        function A = computeLHS(obj)
            npnod = obj.mesh.innerMeshOLD.npnod;
            A = sparse(npnod,npnod);
            for iInt = 1:obj.nInt
                A = A + obj.integrators{iInt}.computeLHS();
            end
        end
        
        function f = integrate(obj,nodalFunc)
            f = cell(1,obj.nInt);
            for iInt = 1:obj.nInt
                f{iInt} = obj.integrators{iInt}.integrate(nodalFunc);
            end
        end
        
        function f = integrateAndSum(obj,nodalFunc)
            f = 0;
            for iInt = 1:obj.nInt
                integrator = obj.integrators{iInt};
                if contains(class(integrator),'Composite')
                    intLocal = integrator.integrateAndSum(nodalFunc);
                    npnod = size(nodalFunc,1);
                    obj.RHScells = zeros(npnod,1);
                    obj.RHSsubcells = intLocal;
                    obj.assembleBoxFaceToGlobal(iInt);
                    int = obj.RHScells;
                else
                    int = integrator.integrate(nodalFunc);
                end
                f = f + int;
            end
        end
        
    end
    
    methods (Access = private)
        
        function createNint(obj,cParams)
            obj.nInt = numel(cParams.compositeParams);
        end
        
        function createIntegrators(obj,cParams)
            obj.createNint(cParams);
            params = cParams.compositeParams;
            for iInt = 1:obj.nInt
                s = params{iInt};
                integrator = Integrator.create(s);
                obj.integrators{end+1} = integrator;
            end
        end
        
        function assembleBoxFaceToGlobal(obj,iInt)
            boxFaceCells = obj.integrators{iInt}.boxFaceToGlobal;
            obj.RHScells(boxFaceCells,:) = obj.RHSsubcells;
        end
        
    end
       
end

