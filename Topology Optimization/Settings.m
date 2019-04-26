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
        constraint_case = 'EQUALITY'
        HJiter0 
        e2 
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
        line_search
        kappaMultiplier
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
            obj.materialInterpolation = method;
            obj.material = materialType;
            obj.initial_case = initial_case;
            obj.cost = cost;
            obj.weights = weights;
            obj.constraint = constraint;
            obj.nconstr = length(constraint);
            obj.optimizer = optimizer;
            obj.kappaMultiplier = kappaMultiplier;
            obj.filter = filterType;
            obj.nsteps = nsteps;
            obj.TOL.rho_plus = TOL.rho_plus;
            obj.TOL.rho_minus = TOL.rho_minus;
            obj.TOL.E_plus = TOL.E_plus;
            obj.TOL.E_minus = TOL.E_minus;
            obj.TOL.nu_plus = TOL.nu_plus;
            obj.TOL.nu_minus = TOL.nu_minus;
            obj.Vfrac_initial = Vfrac_initial;
            obj.optimality_initial  = optimality_initial;
            obj.constr_initial = constr_initial;
            obj.optimality_final = optimality_final;
            obj.constr_final = constr_final;
            
            if exist('line_search','var')
                obj.line_search = line_search;
            else
                if strcmp(obj.optimizer,'PROJECTED GRADIENT')
                    obj.line_search = 'DOUBLING LAST STEP';
                else
                    obj.line_search = 'DIMENSIONALLY CONSISTENT';
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
            
            if ~contains(case_file,'test','IgnoreCase',true)
                fprintf('Loaded %s: \n -Optimizer: %s \n -Cost: ',case_file,obj.optimizer)
                fprintf('%s, ',obj.cost{:})
                fprintf('\n -Constraint: ')
                fprintf('%s, ', obj.constraint{:})
                fprintf('\n -Incremental Steps: %f \n ',obj.nsteps)
            end
            
            if exist('maxiter','var')
                obj.maxiter = maxiter;
                if ~contains(case_file,'test','IgnoreCase',true)
                    fprintf('-Max iters: %f \n ',obj.maxiter)
                end
            end
            
            if exist('Vfrac_final','var')
                obj.Vfrac_final = Vfrac_final;
                if ~contains(case_file,'test','IgnoreCase',true)
                    fprintf('-Volume target: %f \n ',obj.Vfrac_final)
                end
            end
            if exist('Perimeter_target','var')
                obj.Perimeter_target = Perimeter_target;
                if ~contains(case_file,'test','IgnoreCase',true)
                    fprintf('-Perimeter target: %f \n',obj.Perimeter_target)
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