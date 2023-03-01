classdef Sh_volumeColumn < ShapeFunctional
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        mesh
        sectionVariables
        maxVolume
    end
    
    methods (Access = public)
        
        function obj = Sh_volumeColumn(cParams)
            obj.init(cParams)

        end
        function computeFunctionAndGradient(obj)
            obj.computeFunction();
            obj.computeGradient();
        end

        function t = getTitlesToPlot(obj)
            t{1} = 'Volume Column';
        end    

        function v = getVariablesToPlot(obj)
            v{1} = obj.value;%*obj.value0;
        end        

    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.mesh             = cParams.mesh;
            obj.sectionVariables = cParams.sectionVariables;
            obj.maxVolume        = cParams.maxVolume;
        end

    end

    methods (Access = public)

        function computeFunction(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('CONSTANT');
            dV = obj.mesh.computeDvolume(q);
            A = obj.sectionVariables.computeArea();
            V = dV*A;
            Vmax = obj.maxVolume;
            fx = V-Vmax;
            obj.value = fx;
        end

        function computeGradient(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            l = sum(obj.mesh.computeDvolume(q))';              
            nElem = obj.mesh.nelem;
            dA = obj.sectionVariables.computeAreaDerivative();
            nVar = obj.sectionVariables.nDesVarElem;
            dfdx = zeros(1,nVar*nElem+1);
            lVar = repmat(l,nVar,1);
            dfdx(1:nVar*nElem) = dA.*lVar;
            obj.gradient = dfdx;
        end

    end

end