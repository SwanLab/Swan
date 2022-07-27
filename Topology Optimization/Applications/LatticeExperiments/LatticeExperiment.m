classdef LatticeExperiment < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        fileNames
        topOptSet
        topOptProblem
    end
    
    methods (Access = public)
        
        function obj = LatticeExperiment()
            obj.init();
            for icases = 1:numel(obj.fileNames)
                obj.createSettings(icases);
                obj.solveProblem();
                close all
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.fileNames = {   ...
                'CantileverSymmetricFixingDirichletZone'
                %'LshapeCoarseSuperEllipseDesignVariable'
                %'CantileverSymmetricFixingMaxStressZone'
                %'CantileverSymmetricWithoutFixing';
                %'ArchSymSuperEllipseDesignVariable'
   %'LshapeCoarseSuperEllipseDesignVariable'
            %       'ExperimentingPlotSuperEllipseArchDesignVariable'
        %  'ExperimentingPlotSuperEllipseArch';
 %              'ExperimentingPlotRectangle';
        %'ExperimentingPlotSuperEllipse';
            %'CantileverSymmetricMeshSuperEllipsePDEDouble'
          %  'LatticeExperimentInputCantileverSymmetricMeshSuperEllipseP1';
         %   'LatticeExperimentInputCantileverSymmetricMeshRectangleP1'
            %  'ArchTriFineSuperEllipsePDEStressNormP64';
            %    'ArchTriFineRectanglePDEStressNormP64';
               % 'LShapeRectanglePDEStressNormP64';
           % 'LShapeSuperEllipsePDEStressNormP64';
         %  'LatticeExperimentInputCantileverSymmetricMeshSuperEllipsePDE'
           %   'LatticeExperimentInputCantileverSymmetricMeshRectanglePDE';
                %'LatticeExperimentInputCantileverSymmetricMeshSuperEllipsePDE'
              %  'LatticeExperimentInputCantileverSymRectangleDoubleLineSearch'
               % 'LatticeExperimentInputCantileverSymFineRectangle'
         %       'LatticeExperimentInputBulkSymSuperEllipse'
                %'LatticeExperimentInputBulkSymRectangle'
                %'LatticeExperimentInputCantileverTriNewRectangle';
       %         'LatticeExperimentInputCantileverTriFineRectangle'
                %'LatticeExperimentInputLshapeExtendedRectangle';
                %'LatticeExperimentInputLshapeExtendedSuperEllipse';
              %  'LatticeExperimentInputBulkFineRectangle';
              %'LatticeExperimentInputBulkSuperEllipse';

            %   'LatticeExperimentInputLshapeTriFineRectangle';
            %  'LatticeExperimentInputLshapeTriFineSuperEllipse';
            %'LatticeExperimentInputArchTriYSuperEllipse';
            %'LatticeExperimentInputArchTriYRectangle';  
             
        %  'LatticeExperimentInputCantileverTriSuperEllipse';
         %   'LatticeExperimentInputCantileverTriRectangle';
             %  'LatticeExperimentInputCantileverTriFineSuperEllipse';
             %  'LatticeExperimentInputCantileverTriFineRectangle'; 
                         %     'LatticeExperimentInputCantileverTriFineEllipse';

               
               %'LatticeExperimentInputArchTriFineRectangle';
               %'LatticeExperimentInputArchTriFineSuperEllipse';
               %'LatticeExperimentInputLshapeTriFineEllipse';
               %'LatticeExperimentInputArchTriFineEllipse';
               };
            
        end
        
        
        function createSettings(obj,icases)
            s = SettingsTopOptProblem(obj.fileNames{icases});
            obj.topOptSet = s;
        end
        
        function solveProblem(obj)
            obj.topOptProblem = TopOpt_Problem(obj.topOptSet);
            obj.topOptProblem.computeVariables(); 
            obj.topOptProblem.postProcess();
        end
    end
    
end