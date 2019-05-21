classdef ShFunWithElasticPdes < ShapeFunctional
    
    properties (Access = protected)
        interpolation
        physicalProblem
    end
    
    properties (Access = protected)
        regDesignVariable
    end
    
    methods (Access = public)
        
        function computeCostAndGradient(obj)
            obj.nVariables = obj.designVariable.nVariables;
            for i = 1:1
                obj.updateHomogenizedMaterialProperties();
                obj.solvePDEs();
                obj.updateAlpha();
            end
            obj.computeFunctionValue();
            obj.computeGradient();
            obj.normalizeFunctionAndGradient()
        end
        
        function plotAlpha(obj)
            x = obj.physicalProblem.mesh.coord(:,1);
            y = obj.physicalProblem.mesh.coord(:,2);
            conn = obj.physicalProblem.mesh.connec;
            nnode = size(conn,2);
            nelem = size(conn,1);
            xn = zeros(nelem,nnode);
            yn = zeros(nelem,nnode);
            for inode = 1:nnode
                nodes = conn(:,inode);
                xn(:,inode) =  x(nodes);
                yn(:,inode) =  y(nodes);
            end
            xp = mean(xn,2);
            yp = mean(yn,2);
            a1 = obj.designVariable.alpha(1,:);
            a2 = obj.designVariable.alpha(2,:);
            figure();
            quiver(xp,yp,a1,a2) ;
            drawnow
        end
        
        function f = getPhysicalProblems(obj)
            f{1} = obj.physicalProblem;
        end
        
    end
    
    methods (Access = protected)
        
        function createEquilibriumProblem(obj,fileName)
            obj.physicalProblem = FEM.create(fileName);
            obj.physicalProblem.preProcess;
            obj.initPrincipalDirections();
        end
        
        function updateHomogenizedMaterialProperties(obj)
            obj.filterDesignVariable();
            obj.homogenizedVariablesComputer.computeCtensor(obj.regDesignVariable);
        end
        
        function filterDesignVariable(obj)
            nx = length(obj.designVariable.value)/obj.designVariable.nVariables;
            x  = obj.designVariable.value;
            for ivar = 1:obj.nVariables
                i0 = nx*(ivar-1) + 1;
                iF = nx*ivar;
                xs = x(i0:iF);
                xf(:,ivar) = obj.filter.getP0fromP1(xs);
            end
            obj.regDesignVariable = xf;
        end
        
        function computeGradient(obj)
            nelem = obj.physicalProblem.geometry.interpolation.nelem;
            ngaus = obj.physicalProblem.element.quadrature.ngaus;
            nstre = obj.physicalProblem.element.getNstre();
            g = zeros(nelem,ngaus,obj.nVariables);
            for igaus = 1:ngaus
                for istre = 1:nstre
                    for jstre = 1:nstre
                        g(:,igaus,:) = g(:,igaus,:) + obj.updateGradient(igaus,istre,jstre);
                    end
                end
            end
            
            gf = zeros(size(obj.Msmooth,1),obj.nVariables);
            for ivar = 1:obj.nVariables
                gs = squeeze(g(:,:,ivar));
                gf(:,ivar) = obj.filter.getP1fromP0(gs);
            end
            g = obj.Msmooth*gf;
            obj.gradient = g(:);
        end
        
    end
    
    methods (Access = private)
        
        function initPrincipalDirections(obj)
            ndim = obj.physicalProblem.mesh.ndim;
            nelem = obj.physicalProblem.element.nelem;
            alpha0 = zeros(ndim,nelem);
            alpha0(1,:) = 1;
            obj.physicalProblem.variables.principalDirections = alpha0;
            obj.designVariable.alpha = alpha0;
        end
        
        function updateAlpha(obj)
            pD = obj.physicalProblem.variables.principalDirections;
            firstPD = squeeze(pD(:,1,:));
            obj.designVariable.alpha = firstPD;
        end
        
    end
    
    methods (Access = protected, Abstract)
        computeFunctionValue(obj)
        updateGradient(obj)
        solvePDEs(obj)
    end
end