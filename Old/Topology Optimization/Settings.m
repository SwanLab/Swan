classdef Settings %< handle%& matlab.mixin.Copyable
    
    properties %optmizer access
        optimizerSettings
        plotting = true
        showBC = true
        BCscale_factor = 0.10
        rotation_per_it = 0
        printing = true
        printing_physics = false
        monitoring = true
        monitoring_interval = 10
        maxiter = 2000
        constraint_case = {'EQUALITY'}
        HJiter0 
        e2 
        ub = 1;
        lb = 0;
        volumeMicro
        superEllipseRatio
    end
    
    properties %target parameters
        Vfrac_initial
        optimality_initial
        constr_initial
        Vfrac_final
        optimality_final
        constr_final
        epsilon_initial
        epsilon_final
        Perimeter_target
        epsilon_isotropy_initial
        epsilon_isotropy_final
        stressNormExponent_initial
        stressNormExponent_final
    end
    
    properties
        perimeter = struct;
    end
    
    properties %LevelSetCreator
        levelSetDataBase
    end
    
    properties    %topopt access
        ptype
        pdim
        case_file
        filename
        material
        initial_case
        cost
        weights
        constraint
        optimizer
        optimizerUnconstrained        
        line_search_initiator
        incrementFactor
        rate
        filter
        unfitted_mesh_algorithm='DELAUNAY'
        TOL = struct;
        target_parameters = struct;
        nsteps
        micro = struct;
        selectiveC_Cstar
        nconstr
        warningHoleBC
        printIncrementalIter
        printChangingFilter
        printMode = 'DesignAndShapes';
        homegenizedVariablesComputer
        materialInterpolation
        designVariable
        vademecumFileName        
        nelem
        m1
        m2
        alpha0
        rho0
        isDesignVariableFixed
        costDomainNotOptimizable
        constraintDomainNotOptimizable
    end
    
    properties %exploring tests
        shFuncParamsName
    end
    
    methods
        function obj = Settings(case_file)
            run(case_file)
            obj.case_file=case_file;
            obj.filename = filename;
            obj.ptype = ptype;
            
            if exist('method','var')
                obj.materialInterpolation = method;
            end
            
            if exist('materialType','var')
                obj.material = materialType;                
            end
            
            if exist('initial_case','var')
                 obj.initial_case = initial_case;
                 if isequal(initial_case,'full') 
                    obj.levelSetDataBase.type = 'initial_case'; 
                 end
            end            
            
            if exist('isDesignVariableFixed','var')
                femReader = FemInputReader_GiD();
                s = femReader.read(obj.filename);
                coord = s.mesh.coord;                
                obj.isDesignVariableFixed.nodes  = isDesignVariableFixed.nodes(coord);
                obj.isDesignVariableFixed.values = isDesignVariableFixed.values(coord);
            end
            
            if exist('costDomainNotOptimizable','var')
                femReader = FemInputReader_GiD();
                s = femReader.read(obj.filename);
                coord = transpose(s.mesh.computeBaricenter);
                for i = 1:numel(costDomainNotOptimizable)
                    cD = costDomainNotOptimizable{i};
                    if ~isempty(cD)
                        obj.costDomainNotOptimizable{i} = cD(coord);
                    else
                        obj.costDomainNotOptimizable{i} = false(size(coord,1));
                    end
                end
            end
            if exist('constraintDomainNotOptimizable','var')
                femReader = FemInputReader_GiD();
                s = femReader.read(obj.filename);
                coord = transpose(s.mesh.computeBaricenter);
                for i = 1:numel(constraintDomainNotOptimizable)
                    cD = constraintDomainNotOptimizable{i};
                    if ~isempty(cD)
                        obj.constraintDomainNotOptimizable{i} = cD(coord);
                    else
                        obj.constraintDomainNotOptimizable{i} = false(size(coord,1),1);
                    end
                end
            end
            
            obj.cost = cost;
            obj.weights = weights;
            obj.constraint = constraint;
            obj.nconstr = length(constraint);
            obj.optimizer = optimizer;
            if exist('incrementFactor','var')
               obj.incrementFactor = incrementFactor;
            end
            obj.filter = filterType;
            obj.nsteps = nsteps;
            if exist('TOL','var')
                obj.TOL.rho_plus = TOL.rho_plus;
                obj.TOL.rho_minus = TOL.rho_minus;
                obj.TOL.E_plus = TOL.E_plus;
                obj.TOL.E_minus = TOL.E_minus;
                obj.TOL.nu_plus = TOL.nu_plus;
                obj.TOL.nu_minus = TOL.nu_minus;                
            end
            
            obj.Vfrac_initial = Vfrac_initial;
            obj.optimality_initial  = optimality_initial;
            obj.constr_initial = constr_initial;
            obj.optimality_final = optimality_final;
            obj.constr_final = constr_final;

            if exist('line_search_initiator','var')
                obj.line_search_initiator = line_search_initiator;
            else
                if strcmp(obj.optimizer,'PROJECTED GRADIENT')
                    obj.line_search_initiator = 'INCREASING LAST STEP';
                else
                    obj.line_search_initiator = 'STANDARD';
                end
            end
            
            if exist('constraint_case','var')
                obj.constraint_case = constraint_case;
            end
            
            if exist('plotting','var')
                obj.plotting = plotting;
            end
            if exist('showBC','var')
                obj.showBC = showBC;
            end
            if exist('BCscale_factor','var')
                obj.BCscale_factor = BCscale_factor;
            end
            if exist('rotation_per_it','var')
                obj.rotation_per_it = rotation_per_it;
            end
            if exist('printing','var')
                obj.printing = printing;
            end
            if exist('printing_physics','var')
                obj.printing_physics = printing_physics;
            end
            if ~obj.printing && obj.printing_physics
                warning('Physical variables will not be printed.')
            end
            if exist('monitoring','var')
                obj.monitoring = monitoring;
            end
            if exist('monitoring_interval','var')
                obj.monitoring_interval = monitoring_interval;
            end
            
%             if ~contains(case_file,'test','IgnoreCase',true)
%                 fprintf('Loaded %s: \n -Optimizer: %s \n -Cost: ',case_file,obj.optimizer)
%                 fprintf('%s, ',obj.cost{:})
%                 fprintf('\n -Constraint: ')
%                 fprintf('%s, ', obj.constraint{:})
%                 fprintf('\n -Incremental Steps: %f \n ',obj.nsteps)
%             end
            
            if exist('maxiter','var')
                obj.maxiter = maxiter;
                if ~contains(case_file,'test','IgnoreCase',true)
%                     fprintf('-Max iters: %f \n ',obj.maxiter)
                end
            end
            
            if exist('Vfrac_final','var')
                obj.Vfrac_final = Vfrac_final;
                if ~contains(case_file,'test','IgnoreCase',true)
%                     fprintf('-Volume target: %f \n ',obj.Vfrac_final)
                end
            end
            if exist('Perimeter_target','var')
                obj.Perimeter_target = Perimeter_target;
                if ~contains(case_file,'test','IgnoreCase',true)
%                     fprintf('-Perimeter target: %f \n',obj.Perimeter_target)
                end
            end
            if exist('epsilon_initial','var')
                obj.epsilon_initial = epsilon_initial;
                obj.epsilon_final = epsilon_final;
            end
            if exist('HJiter0','var')
                obj.HJiter0 = HJiter0;
            end
            if exist('e2','var')
                obj.e2 = e2;
            else
                obj.e2 = 1;
            end

            
            if exist('micro','var')
                obj.micro.alpha = micro.alpha;
                obj.micro.beta = micro.beta;
            end
            if exist('epsilon_isotropy_initial','var')
                obj.epsilon_isotropy_initial = epsilon_isotropy_initial;
                obj.epsilon_isotropy_final = epsilon_isotropy_final;
            end
            if exist('selectiveC_Cstar','var')
                obj.selectiveC_Cstar = selectiveC_Cstar;
            end
            
            if exist('widthH','var')
                obj.levelSetDataBase.widthH = widthH;
            end
            
            if exist('widthV','var')
                obj.levelSetDataBase.widthV = widthV;
            end
            
            if exist('widthSquare','var')
                obj.levelSetDataBase.widthSquare = widthSquare;
            end
            
            if exist('N_holes','var')
                obj.levelSetDataBase.nHoles = N_holes;
            end
            if exist('R_holes','var')
                obj.levelSetDataBase.rHoles = R_holes;
            end
            if exist('phase_holes','var')
                obj.levelSetDataBase.phaseHoles = phase_holes;
            end  
            
            if exist('warningHoleBC','var')
                obj.levelSetDataBase.warningHoleBC = warningHoleBC;
            end  
            
            if exist('fracRadius','var')
                obj.levelSetDataBase.fracRadius = fracRadius;
            end
            
            if exist('levFib','var')
                obj.levelSetDataBase.levFib = levFib;
            end  
            
            if exist('yn','var')
                obj.levelSetDataBase.yn = yn;
            end   
            
            if exist('levelFibers','var')
               obj.levelSetDataBase.levelFibers = levelFibers;               
            end
            
            if exist('volumeFibers','var')
               obj.levelSetDataBase.volumeFibers = volumeFibers;               
            end
            
            if exist('shFuncParamsName','var')
               obj.shFuncParamsName = shFuncParamsName;               
            end
            
            
            
            if exist('designVariable','var')
               obj.designVariable = designVariable;               
            end
            
            if exist('homegenizedVariablesComputer','var')
               obj.homegenizedVariablesComputer = homegenizedVariablesComputer;               
            else
               obj.homegenizedVariablesComputer = 'ByInterpolation'; 
               
            end

            if exist('vademecumFileName','var')            
                obj.vademecumFileName = vademecumFileName;
            end
            
            
            if  ~(contains(filename,'test','IgnoreCase',true) || contains(filename,'RVE') || obj.hasToAddSpaceBecauseOfIncremental())
                fprintf('\n')
            end
            
            if ~exist('optimizerUnconstrained','var')
                switch obj.optimizer
                    case {'SLERP','HAMILTON-JACOBI','PROJECTED GRADIENT'}
                        obj.optimizerUnconstrained = obj.optimizer;
                        obj.optimizer = 'AlternatingPrimalDual';
                end
            else                
                obj.optimizerUnconstrained = optimizerUnconstrained;
            end
            
            if exist('rate','var')
                obj.rate = rate;
            else 
                obj.rate = 0.5;
            end
            
            if exist('ub','var')
                obj.ub = ub;
            else
                obj.ub = 1;
            end
            
            if exist('lb','var')
                obj.lb = lb;
            else
                obj.lb = 0;
            end   
            
            if exist('superEllipseRatio','var')
               obj.levelSetDataBase.vigdergauzDataBase.superEllipseRatio = superEllipseRatio;
            end
            
            if exist('volumeMicro','var')
               obj.levelSetDataBase.vigdergauzDataBase.volumeMicro = volumeMicro;
            end            
            
            if exist('vigdergauzType','var')
               obj.levelSetDataBase.vigdergauzDataBase.type =  vigdergauzType;
            end
            
            if exist('vigdergauzStrainMacro','var')
               obj.levelSetDataBase.vigdergauzDataBase.strain =  vigdergauzStrainMacro;
               if exist('TOL','var')
                   obj.levelSetDataBase.vigdergauzDataBase.E1 = TOL.E_plus;
                   obj.levelSetDataBase.vigdergauzDataBase.E0 = TOL.E_minus;
                   obj.levelSetDataBase.vigdergauzDataBase.nu1 = TOL.nu_plus;
                   obj.levelSetDataBase.vigdergauzDataBase.nu0 = TOL.nu_minus;
               end
            end   
            
            if exist('rho0','var')
                obj.rho0 = rho0;
            end
            
            if exist('m1','var')
               obj.m1 = m1;
            end
            
            if exist('m2','var')
               obj.m2 = m2;
            end  
            
            if exist('alpha0','var')
               obj.alpha0 = alpha0;
            end
            
            if exist('stressNormExponent_initial','var')
                obj.stressNormExponent_initial = stressNormExponent_initial;
            end            
            
            if exist('stressNormExponent_final','var')
                obj.stressNormExponent_final = stressNormExponent_final;
            end            
            
        end
        
       
    end
    
    
    methods (Access = private)
        
        function itHas = hasToAddSpaceBecauseOfIncremental(obj)
            itHas = true;
            if ~isempty(obj.printIncrementalIter)
                itHas = obj.printIncrementalIter;
            end
        end
        
    end
end