classdef Sh_volumeColumn < ShapeFunctional
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        mesh
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
            obj.designVariable = cParams.designVariable;
            obj.mesh           = cParams.mesh;
        end

    end

    methods (Access = public)

        function computeFunction(obj)
            V = obj.designVariable.computeVolum();
            fx = V-1;
           % fx = sqrt(V)-1;
            obj.value = fx;
        end

        function computeGradient(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            l = sum(obj.mesh.computeDvolume(q));             
            %dfdx(1,:) = 1./1.*dfdx(1,:); 
            %dfdx(2,:) = 1./1.*dfdx(2,:); 
            nElem = obj.designVariable.mesh.nelem;
            dfdx = zeros(1,nElem+1);
            dfdx(1,1:nElem)= l.';%*ones(1,nElem);% 1/(nElem+1).*ones(1,nElem);
            obj.gradient = dfdx;
        end

    end

end