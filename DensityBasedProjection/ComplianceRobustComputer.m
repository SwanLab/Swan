classdef ComplianceRobustComputer < handle
    properties (Access = public)
        projectedField        
    end
    properties (Access = private)
        iterations
        filteredField
        field
        mesh
        structure
        projectParameters
        filterParameters
        solverParameters
        cost
        volumen
    end
    methods (Access = public)
        function obj = ComplianceRobustComputer(cParams)
            obj.inputData(cParams);
        end
        function compute(obj)
            obj.computeInitialParameters();
            obj.projectField();
            obj.computeCost();
            obj.optimize();
        end
    end
    methods (Access = private)
        function inputData(obj,cParams)
            obj.iterations = cParams.iterations;
        end 
        function computeInitialParameters(obj)
            obj.computeMeshParameters();
            obj.computeStructureParameters();
            obj.computeSolverParameters();
            obj.computeProjectorParameters
            obj.computeFilterParameters();
            obj.computeInitialFields();
        end
        function projectField(obj)
            % Project the initial field
            s.beta = obj.projectParameters.beta;
            s.filteredField =obj.filteredField;
            s.eta =obj.projectParameters.eta.E;
            E = FieldProjector(s);
            E.compute();
            obj.projectedField.E = E.projectedField;

            s.eta =obj.projectParameters.eta.I;
            I = FieldProjector(s);
            I.compute();
            obj.projectedField.I = I.projectedField;

            s.eta =obj.projectParameters.eta.D;
            D = FieldProjector(s);
            D.compute();
            obj.projectedField.D = D.projectedField;
            obj.volumen.volfracD = obj.filterParameters.volumenFraction*sum(obj.projectedField.D(:))/sum(obj.projectedField.I(:));
        end 
        function computeCost(obj)          
            %Get intial cost
            s.mesh = obj.mesh; 
            s.structure = obj.structure;
            s.projectedField = obj.projectedField.E;

            B = FEMcomputer(s);
            B.compute();
            obj.cost.E = B.cost;
        end 
        function optimize(obj)
            %Optimización
            s.mesh = obj.mesh;
            s.structure=obj.structure;
            s.structure.elementType = 'Square';
            s.projector = obj.projectParameters;
            s.filterParameters = obj.filterParameters;
            s.cost = obj.cost;
            s.cost.initial = obj.cost.E;
            s.solverParameters =obj.solverParameters;
            s.iterations = obj.iterations;
            s.field = obj.field;
            s.filteredField = obj.filteredField;
            s.projectedField = obj.projectedField;
            s.volumen = obj.volumen;
            B = Optimizer(s);
            B.compute();
            obj.projectedField = B.projectedField;
        end 
        function computeMeshParameters(obj)
            %Mesh parameters
            %Malla estandar 400x80
            obj.mesh.elementNumberX    = 160;
            obj.mesh.elementNumberY    = 80;
            obj.mesh.neumanCondition       = -1e-3;
            obj.mesh.output    = 2;

            s.elementNumberX =  obj.mesh.elementNumberX;
            s.elementNumberY =  obj.mesh.elementNumberY;
            B = GeometryComputer(s);
            B.compute();
            obj.mesh.degress.all   = B.degress.all;
            obj.mesh.degress.free  = B.degress.free;
            obj.mesh.degress.fixed = B.degress.fixed;
            obj.mesh.conectivityMatrixMat = B.conectivityMatrixMat;
        end
        function computeStructureParameters(obj)
             obj.structure.t       = 1;
             %Parámetros del material:
            obj.structure.elasticModuleNeutral      = 1;
            obj.structure.elasticModuleMinimun    = 1e-6;
            obj.structure.poissonCoefficient      = 0.3;
            obj.structure.penalization   = 3;

            %Compute Elemental stifness Matrix
            s.elementType = 'square';
            s.t = obj.structure.t;
            s.poissonCoefficient =obj.structure. poissonCoefficient;
            B = ElementalStiffnessMatricesComputer(s);
            B.compute();
            obj.structure.elementalStiffnessMatrix = B.elementalStiffnessMatrix;
        end 
        function computeSolverParameters(obj)
            obj.cost.change  = 1;
            obj.solverParameters.minDensity    = zeros(obj.mesh.elementNumberY,obj.mesh.elementNumberX);
            obj.solverParameters.maxDensity    =  ones(obj.mesh.elementNumberY,obj.mesh.elementNumberX);
            obj.solverParameters.minDensity =  obj.solverParameters.minDensity(:);
            obj.solverParameters.maxDensity =  obj.solverParameters.maxDensity(:);
            
            obj.solverParameters.numberRestriction       = 4;
            obj.solverParameters.variableNumber       = obj.mesh.elementNumberX*obj.mesh.elementNumberY;
            obj.solverParameters.low     = obj.solverParameters.minDensity;
            obj.solverParameters.upp     = obj.solverParameters.maxDensity;
            obj.solverParameters.a0      = 1;
            obj.solverParameters.mmaParameter.a       = [1 1 1 0]';
            obj.solverParameters.mmaParameter.e       = 1000*ones(obj.solverParameters.numberRestriction,1);
            obj.solverParameters.mmaParameter.d       = 0*ones(obj.solverParameters.numberRestriction,1);
            obj.solverParameters.xold1   = obj.field(:);
            obj.solverParameters.xold2   = obj.field(:);
        end 
        function computeFilterParameters(obj)
            %Filter parameters:
            obj.filterParameters.volumenFraction =   0.5;
            obj.filterParameters.minimunInfluenceRadios    =   1.5;
            
            s.elementNumberX =  obj.mesh.elementNumberX;
            s.elementNumberY =  obj.mesh.elementNumberY;
            s.minimunInfluenceRadios = obj.filterParameters.minimunInfluenceRadios;
            B = weightFilterComputer(s);
            B.compute();
            obj.filterParameters.H = B.H;
            obj.filterParameters.Hs = B.Hs;
        end
        function computeProjectorParameters(obj)
            eta      = 0.25;
            obj.projectParameters.eta.E     = 1-eta;
            obj.projectParameters.eta.I     = 0.5;
            obj.projectParameters.eta.D     = eta;
            obj.projectParameters.beta     = 1;
        end
        function computeInitialFields(obj)
            %Define initial fields
            obj.field        = obj.filterParameters.volumenFraction*ones(obj.mesh.elementNumberY,obj.mesh.elementNumberX);
            obj.filteredField = obj.field;
            %Define density constrains
            obj.field([],[]) = 1;
        end
        
    end
end