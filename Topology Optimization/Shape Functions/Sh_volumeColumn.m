classdef Sh_volumeColumn < ShapeFunctional
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        mesh
        sectionVariables
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
            %obj.designVariable = cParams.designVariable;
            obj.mesh           = cParams.mesh;
            obj.sectionVariables = cParams.sectionVariables;
        end

    end

    methods (Access = public)

        function computeFunction(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('CONSTANT');
            dV = obj.mesh.computeDvolume(q);
            A = obj.sectionVariables.computeArea();
            V = dV*A;
%             V = obj.designVariable.computeVolum();
            fx = V-1;
%           fx = sqrt(V)-1;
            obj.value = fx;
        end

        function computeGradient(obj)
            %R = obj.designVariable.getColumnRadius();
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            l = sum(obj.mesh.computeDvolume(q));             
            %dfdx(1,:) = 1./1.*dfdx(1,:); 
            %dfdx(2,:) = 1./1.*dfdx(2,:); 
            nElem = obj.mesh.nelem;
            dfdx = zeros(1,2*nElem+1);
%             dfdx = zeros(1,nElem+1);
            %dfdx(1,1:nElem)= 2*pi*l'.*R; %*ones(1,nElem);% 1/(nElem+1).*ones(1,nElem);
            dA = obj.sectionVariables.computeAreaDerivative();
            %dfdx(1,1:nElem)= dA.*l';
            dfdx(1,1:2*nElem)= dA.*[l,l]';
            obj.gradient = dfdx;
        end

    end

end