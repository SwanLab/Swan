classdef ComplianceRobustComputer < handle
    properties (Access = public)
        projectedField
        iterations
    end
    properties (Access = private)
    end
    methods (Access = public)
        function obj = ComplianceRobustComputer(cParams)
            obj.inputData(cParams);
        end
        function compute(obj)
            obj.runCode();
        end
    end
    methods (Access = private)
        function inputData(obj,cParams)
            obj.iterations = cParams.iterations;
        end 
        
        function runCode(obj)

            %% Parameters

            %Malla (elementos cuadrados) y dimensiones de la estructura:
            %Malla estandar 400x80
            elementNumberX    = 160;
            elementNumberY    = 80;
            t       = 1;

            %Parámetros del material:
            elasticModuleNeutral      = 1;
            elasticModuleMinimun    = 1e-6;
            poissonCoefficient      = 0.3;
            penalization   = 3;

            %Entradas del problema:
            volumenFraction =   0.5;
            minimunInfluenceRadios    =   1.5;
            neumanCondition       = -1e-3;
            output    = 2;

            %Define initial fields
            field        = volumenFraction*ones(elementNumberY,elementNumberX);
            filteredField = field;

            %Define density constrains
            field([],[]) = 1;
            minDensity    = zeros(elementNumberY,elementNumberX);
            maxDensity    =  ones(elementNumberY,elementNumberX);
            minDensity = minDensity(:);
            maxDensity = maxDensity(:);

            %MMA parameters
            costChange  = 1;
            numberRestriction       = 4;
            variableNumber       = elementNumberX*elementNumberY;
            low     = minDensity;
            upp     = maxDensity;
            a0      = 1;
            mmaParameter.a       = [1 1 1 0]';
            mmaParameter.e       = 1000*ones(numberRestriction,1);
            mmaParameter.d       = 0*ones(numberRestriction,1);
            xold1   = field(:);
            xold2   = field(:);

            %Compute Elemental stifness Matrix
            s.elementType = 'square';
            s.t = t;
            s.poissonCoefficient = poissonCoefficient;
            B = ElementalStiffnessMatricesComputer(s);
            B.compute();
            elementalStiffnessMatrix = B.elementalStiffnessMatrix;

            %Mesh parameters
            s.elementNumberX =  elementNumberX;
            s.elementNumberY =  elementNumberY;
            B = GeometryComputer(s);
            B.compute();
            allDegress   = B.degress.all;
            freeDegress  = B.degress.free;
            fixedDegress = B.degress.fixed;
            conectivityMatrixMat = B.conectivityMatrixMat;


            %Filter parameters:
            s.elementNumberX =  elementNumberX;
            s.elementNumberY =  elementNumberY;
            s.minimunInfluenceRadios = minimunInfluenceRadios;
            B = weightFilterComputer(s);
            B.compute();
            H = B.H;
            Hs = B.Hs;

            eta      = 0.25;
            etaE     = 1-eta;
            etaI     = 0.5;
            etaD     = eta;
            beta     = 1;


            % Project the initial field
            s.beta = beta;
            s.eta =etaE;
            s.filteredField =filteredField;
            E = FieldProjector(s);
            E.compute();
            projectedFieldE = E.projectedField;

            s.beta = beta;
            s.eta =etaI;
            s.filteredField =filteredField;
            I = FieldProjector(s);
            I.compute();
            projectedFieldI = I.projectedField;

            s.beta = beta;
            s.eta =etaD;
            s.filteredField =filteredField;
            D = FieldProjector(s);
            D.compute();
            projectedFieldD = D.projectedField;

            volfracD = volumenFraction*sum(projectedFieldD(:))/sum(projectedFieldI(:));

            %Get intial cost
            s.mesh.elementNumberX = elementNumberX;
            s.mesh.elementNumberY = elementNumberY;
            s.mesh.neumanCondition = neumanCondition;
            s.mesh.output = output;
            s.mesh.freeDegress = freeDegress;
            s.mesh.conectivityMatrixMat = conectivityMatrixMat;

            s.structure.elementalStiffnessMatrix = elementalStiffnessMatrix;
            s.structure.t=t;
            s.structure.penalization=penalization;
            s.structure.poissonCoefficient=poissonCoefficient;
            s.structure.elasticModuleMinimun=elasticModuleMinimun;
            s.structure.elasticModuleNeutral=elasticModuleNeutral;
            s.projectedField = projectedFieldE;

            B = FEMcomputer(s);
            B.compute();
            cost.E = B.cost;


            %Optimización
            clear s
            s.field = field;
            s.filteredField =filteredField;
            s.mesh.elementNumberX = elementNumberX;
            s.mesh.elementNumberY = elementNumberY;
            s.mesh.neumanCondition = neumanCondition;
            s.mesh.output =output;
            s.mesh.freeDegress = freeDegress;
            s.mesh.conectivityMatrixMat = conectivityMatrixMat;

            s.structure.elementType = 'Square';
            s.structure.t=t;
            s.structure.penalization=penalization;
            s.structure.poissonCoefficient=poissonCoefficient;
            s.structure.elasticModuleMinimun=elasticModuleMinimun;
            s.structure.elasticModuleNeutral=elasticModuleNeutral;
            s.structure.elementalStiffnessMatrix=elementalStiffnessMatrix;

            s.projector.eta.E = etaE;
            s.projector.eta.I = etaI;
            s.projector.eta.D = etaD;
            s.projector.beta = beta;

            s.filter.H = H;
            s.filter.Hs = Hs;

            s.cost.initial = cost.E;
            s.costChange = costChange;

            s.solver.minDensity =minDensity;
            s.solver.maxDensity = maxDensity;
            s.solver.mmaParameter =mmaParameter;
            s.solver.numberRestriction = numberRestriction;
            s.solver.variableNumber = variableNumber;
            s.solver.xold1 = xold1;
            s.solver.xold2 = xold2;
            s.solver.low     = low;
            s.solver.upp     = upp;
            s.solver.a0      = a0;



            s.iterations = obj.iterations;

            s.field = field;
            s.filteredField = filteredField;

            s.projectedField.E = projectedFieldE;
            s.projectedField.I = projectedFieldI;
            s.projectedField.D = projectedFieldD;

            s.problemParameters.volfracD = volfracD;


            
            B = Optimizer(s);
            B.compute();
            obj.projectedField = B.projectedField;
        end
    end
end