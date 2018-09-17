classdef testHomogenization < handle
    
    %UNTITLED8 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ChChecked
        ChToCheck
        fileName
        Interpolation
    end
    
    
    methods
        
        function obj = testHomogenization()
            obj.printHeader()
            obj.computeNumericallyChForLaminate();
            obj.computeSequentialLaminateTensor()
            obj.checkTestPassed()
            obj.printTail()
        end
        
        function computeSequentialLaminateTensor(obj)
            %SeqLam = SequentialLaminate();
            %obj.ChToCheck = SeqLam.Ch;
            Fraction = 0.8;
            
            matPropsStiff = obj.Interpolation.computeMatProp(1);
            MpStiff = Material_Elastic_ISO_2D(1);
            MpStiff = MpStiff.setProps(matPropsStiff);
            Cstiff = MpStiff.C;
            
            
            matPropsWeak.mu = 1e-3*matPropsStiff.mu;
            matPropsWeak.kappa = 1e-3*matPropsStiff.kappa;
            MpWeak = Material_Elastic_ISO_2D(1);
            MpWeak =  MpWeak.setProps(matPropsWeak);            
            Cweak = MpWeak.C;
            
            Chomog = (Fraction)*Cstiff + (1 - Fraction)*Cweak;
            
            d1 = 1; d2 = 0; d3 = 0;
            direction = [d1 d2 d3]; 
            
            SeqLam = SequentialLaminateHomogenizer(matPropsStiff.mu,matPropsStiff.kappa, ...
                                                  matPropsWeak.mu,matPropsWeak.kappa,direction);

            obj.ChToCheck = Chomog;
        end
        
        
        function computeNumericallyChForLaminate(obj)
            MicroTest = {'test_microHorizontal'};
            obj.fileName = MicroTest{1};
            %file_name_in = strcat('./Input/',MicroTest);
            path = '.';
            load_file = strcat(path,'/tests/',obj.fileName);            
            load(load_file);
            file_name_in = strcat('./Input/',MicroTest{1});
            settings = Settings(file_name_in);
            settings.pdim = '2D';
            obj.Interpolation = Material_Interpolation.create(settings.TOL,settings.material,settings.method,settings.pdim);
            Filter = Filter_P1_LevelSet_2D(settings.filename,settings.ptype);
            
            
            MicroProblem = Elastic_Problem_Micro(settings.filename);
            MicroProblem.preProcess();

            
            design_variable_initializer = DesignVaribleInitializer.create(settings,MicroProblem.mesh,MicroProblem.mesh.mean_cell_size);
            design_variable_initializer.compute_initial_design();
            
            %plot_nodal_field(design_variable_initializer.x,MicroProblem.mesh.coord)
            Filter.preProcess();
            rho= Filter.getP0fromP1(design_variable_initializer.x);
            matProps= obj.Interpolation.computeMatProp(rho);
            MicroProblem.setMatProps(matProps);
            
            Vol = ShFunc_Volume(settings);
            Vol.filter.preProcess();
            Vol.computeCostAndGradient(design_variable_initializer.x)
            
            MicroProblem.computeChomog();
            obj.ChChecked = MicroProblem.variables.Chomog;
        end
        
        
        function hasPassed = hasPassed(obj)
            hasPassed = norm(obj.ChChecked - obj.ChToCheck) < 1e-6;
        end
        
        function checkTestPassed(obj)
            if obj.hasPassed()
                cprintf('green',strcat(obj.fileName,' PASSED\n'));
            else
                cprintf('err',strcat(obj.fileName,' FAILED\n'));
            end
        end
        
        
    end
    
    methods (Static)
        
        
 
        function printHeader()
            fprintf('Running TopOpt tests...\n')
        end
        
        function printTail()
            fprintf('\nTopOpt tests completed.\n')
            fprintf('\n-------------------------------------------\n\n')
        end
        

    end
end