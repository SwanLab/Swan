classdef ShFunWithElasticPdes < ShapeFunctional
    
    properties (Access = protected)
        interpolation
        physicalProblem
    end
    
    properties (Access = private)
        regDesignVariable
    end
    
    methods (Access = public)
        
        function computeCostAndGradient(obj)
            obj.nVariables = obj.designVariable.nVariables;
            for i = 1:1
                obj.updateHomogenizedMaterialProperties();
                obj.solvePDEs();
                %obj.initPrincipalDirections();                
                obj.updateAlpha();
                obj.plotAlpha();
            end
            obj.computeFunctionValue();
            obj.computeGradient();
            obj.normalizeFunctionAndGradient()
        end
        
        function plotAlpha(obj)
            x = obj.physicalProblem.mesh.coord(:,1);
            y = obj.physicalProblem.mesh.coord(:,2);
            conn = obj.physicalProblem.mesh.connec;
            for inode = 1:3
                nodes = conn(:,inode);
                xp(:,inode) =  x(nodes);
                yp(:,inode) =  y(nodes);
            end
            xp = mean(xp');
            yp = mean(yp');
            a1 = obj.designVariable.alpha(1,:);
            a2 = obj.designVariable.alpha(2,:);
            figure(200);
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
            ndim = 2;
            nelem = obj.physicalProblem.element.nelem;
            alpha0 = zeros(ndim,nelem);
            alpha0(1,:) = 1;
            obj.physicalProblem.variables.principalDirections = alpha0;
            obj.designVariable.alpha = alpha0;
            obj.designVariable.beta = alpha0;
        end
        
        function updateAlpha(obj)
            stress = squeeze(obj.physicalProblem.variables.avarageStress);
            sx = stress(1,:);
            sy = stress(2,:);
            sxy = stress(3,:);
            sn = sqrt((sx - sy).^2 + 4*sxy.^2);
            b(1,:) = (sx - sy)./sn;
            b(2,:) = 2*sxy./sn;
            
            pD = obj.physicalProblem.variables.principalDirections;
            pS = obj.physicalProblem.variables.principalStress;
            [~,indM] = min(abs(pS));
            for i = 1:2
                dirD = squeeze(pD(i,:,:));
                for j = 1:2
                    ind = indM == j;
                    dir(i,ind) = dirD(j,ind);
                end
            end
            obj.designVariable.alpha = squeeze(pD(:,1,:));%dir;
            obj.designVariable.beta = b;
        end
        
    end
    
    methods (Access = protected, Abstract)
        computeFunctionValue(obj)
        updateGradient(obj)
        solvePDEs(obj)
    end
end