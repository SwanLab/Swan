classdef CoarsePlotSolution
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Property1
    end

    methods
        function obj = CoarsePlotSolution(x, mesh, bcApplier,outputFileName)
            row = 1;
            col = 1;
            iter = 1;
            flag = 0;

            if ~isempty(bcApplier)
                x = bcApplier.reducedToFullVectorDirichlet(x);
            end
            if nargin <7
                flag =0;
            end
            %             xFull = bc.reducedToFullVector(x);
            if size(x,2)==1
                s.fValues = reshape(x,2,[])';
            else
                s.fValues = x;
            end
            %

            s.mesh = mesh;
            s.fValues(:,end+1) = 0;
            s.ndimf = 2;
            s.order = 'P1';
            xF = LagrangianFunction(s);
            %             xF.plot();
            if flag == 0
                xF.print(['domain',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
            elseif flag == 1
                xF.print(['DomainResidual',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
            elseif flag == 2
                xF.print(['Residual',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
            elseif flag == 3
                xF.print(['domainFine',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
            elseif flag == 4
                xF.print(['domainNeuman',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
            end
            
            fileName = ['domain',num2str(row),num2str(col),'_',num2str(iter)];

            s = dir(pwd);
            s = struct2table(s);
            idx = startsWith(s.name, fileName);
            s = s(idx,:);
            oldFileName = s.name;
            newFileName = replace(oldFileName, fileName, outputFileName);
            fclose('all');
            movefile(oldFileName{1}, newFileName{1})

            
        end

        
    end
end