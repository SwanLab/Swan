classdef BistableSpecimenComputer < handle

    properties (GetAccess = public, SetAccess = private)
        geometry
        mesh
        results
    end

    properties (Access = private)
        name
        refinement
        conditionsType
        matProp
        strains
        tolerance
        maxIter
        quadOrder
        outDir

        model
        conditions
        functional
        energyEvaluator
        displacements
    end

    methods (Access = public)

        function obj = BistableSpecimenComputer(cParams)
            obj.init(cParams);
            obj.createMesh();
            obj.createConditions();
        end

        function compute(obj)
            obj.createFunctional();
            output = obj.solve();
            obj.displacements = output.displacement.function;
            obj.printDisplacements();
            obj.results = obj.computeResults(output);
        end

        function plotSetup(obj)
            f = figure('Visible','off','Position',[0 0 1000 900]);
            obj.mesh.plot();
            axis equal; hold on
            x     = obj.conditions.nodeCoordinates;
            nodes = obj.computeRestrictedNodes();
            plot(x(nodes,1),x(nodes,2),'o','MarkerSize',3,'MarkerFaceColor','b','MarkerEdgeColor','b');
            title(sprintf('%s: %d nodos, %d elementos, %d nodos restringidos', ...
                obj.name,size(x,1),obj.mesh.nelem,numel(nodes)),'Interpreter','none');
            xlabel('x [mm]'); ylabel('y [mm]');
            exportgraphics(f,fullfile(obj.outDir,[obj.name,'_setup.png']),'Resolution',130);
            close(f);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.name           = cParams.name;
            obj.geometry       = cParams.geometry;
            obj.refinement     = cParams.refinement;
            obj.conditionsType = cParams.conditionsType;
            obj.matProp        = cParams.matProp;
            obj.strains        = cParams.strains;
            obj.tolerance      = cParams.tolerance;
            obj.maxIter        = cParams.maxIter;
            obj.quadOrder      = 3;
            obj.outDir         = cParams.outDir;
            obj.model          = SpecimenModel.create(cParams);
        end

        function createMesh(obj)
            obj.mesh = obj.model.createMesh(obj.refinement,obj.outDir,obj.name);
        end

        function createConditions(obj)
            s.type     = obj.conditionsType;
            s.values   = obj.strains;
            s.geometry = obj.geometry;
            s.order    = obj.model.interpolationOrder();
            obj.conditions = SpecimenConditionsCreator(obj.mesh,s);
        end

        function createFunctional(obj)
            p = obj.model.createMaterialProperties(obj.matProp,obj.mesh);
            s.quadOrder   = obj.quadOrder;
            s.matProp     = p;
            s.mesh        = obj.mesh;
            s.testSpace.u = LagrangianFunction.create(obj.mesh,2,obj.model.interpolationOrder());
            s = obj.model.addContact(s);
            obj.functional      = ElasticityFunctional(s);
            obj.energyEvaluator = obj.model.createEnergyEvaluator(p,obj.mesh);
        end

        function output = solve(obj)
            s.mesh               = obj.mesh;
            s.boundaryConditions = obj.conditions;
            s.functional         = obj.functional;
            s.material           = obj.createSolidMaterial();
            s.monitoring         = obj.createMonitoringSettings();
            s.tolerance          = obj.tolerance;
            s.maxIter            = obj.maxIter;
            s.stopWhenNotConverged = true;
            s.usePredictor       = true;
            s.lineSearch         = true;
            s.residualTolerance  = 1e-6;
            s.viscosity          = 1e-3;
            s.checkOverlaps      = obj.model.checksOverlaps();
            s.interpolationOrder = obj.model.interpolationOrder();
            s.maxBisections      = 4;
            output = HyperelasticityComputer(s).compute();
        end

        function m = createSolidMaterial(obj)
            N      = obj.mesh.ndim;
            mu     = ConstantFunction.create(obj.matProp.mu,obj.mesh);
            lambda = ConstantFunction.create(obj.matProp.lambda,obj.mesh);
            I      = ConstantFunction.create(eye4D(N),obj.mesh);
            IxI    = ConstantFunction.create(kronEye(N),obj.mesh);
            m.tensor      = 2*Expand(mu,4).*I + Expand(lambda,4).*IxI;
            m.prop.mu     = mu;
            m.prop.lambda = lambda;
        end

        function s = createMonitoringSettings(obj)
            s.set         = false;
            s.printInfo   = true;
            s.printFile   = false;
            s.fileNameOut = obj.name;
        end

        function r = computeResults(obj,output)
            n = numel(obj.displacements);
            r.name     = obj.name;
            r.type     = obj.conditionsType;
            r.strain   = obj.strains(1:n)';
            r.iter     = output.iter(:);
            r.overlaps = obj.readOverlaps(output,n);
            [r.energy,r.force,r.residual,r.minJ,r.maxStretch] = deal(zeros(n,1));
            for i = 1:n
                [r.energy(i),r.force(i),r.residual(i),r.minJ(i),r.maxStretch(i)] = ...
                    obj.computeStepResults(obj.displacements{i});
            end
            r.nNodes = obj.mesh.nnodes;
            r.nElem  = obj.mesh.nelem;
        end

        function [energy,force,residual,minJ,maxStretch] = computeStepResults(obj,u)
            fInt     = obj.energyEvaluator.computeGradient(u,obj.quadOrder);
            energy   = obj.energyEvaluator.computeCost(u,obj.quadOrder);
            force    = obj.conditions.computeConjugateForce(fInt);
            residual = norm(fInt(obj.computeFreeDofs(u)));
            [minJ,maxStretch] = obj.computeDeformationBounds(u);
        end

        function free = computeFreeDofs(obj,u)
            free = setdiff((1:numel(u.fValues))',obj.conditions.computeRestrictedDofs());
        end

        function nodes = computeRestrictedNodes(obj)
            nodes = unique(ceil(obj.conditions.computeRestrictedDofs()/2));
        end

        function printDisplacements(obj)
            for i = 1:numel(obj.displacements)
                s.mesh     = obj.mesh;
                s.type     = 'Paraview';
                s.filename = fullfile(obj.outDir,sprintf('%s_u_step%03d',obj.name,i));
                [s.fun,s.funNames] = obj.model.functionsToPrint(obj.displacements{i});
                FunctionPrinter.create(s).print();
            end
        end

    end

    methods (Static, Access = private)

        function overlaps = readOverlaps(output,n)
            overlaps = zeros(n,1);
            if isfield(output,'overlaps')
                overlaps = output.overlaps(:);
            end
        end

        function [minJ,maxStretch] = computeDeformationBounds(u)
            F      = Identity(u) + Grad(u)';
            jac    = Det(F);
            trC    = DDP(F,F);
            lambda = ((trC + (trC.^2 - 4.*jac.^2).^0.5)./2).^0.5;
            centroid   = [1/3;1/3];
            minJ       = min(jac.evaluate(centroid),[],'all');
            maxStretch = max(real(lambda.evaluate(centroid)),[],'all');
        end

    end

end
