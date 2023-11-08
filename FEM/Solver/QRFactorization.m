classdef QRFactorization < Solver

    methods(Static)

        function [Q, R] = solve(LHS,RHS,mesh,bc)

            [m, n] = size(LHS);
            Q = zeros(m, n);
            R = zeros(n, n);
            x = zeros(m, 1);
   
            for j = 1:n
                v = LHS(:, j);
                for i = 1:j-1
                    R(i, j) = Q(:, i)' * LHS(:, j);
                    v = v - R(i, j) * Q(:, i);
                end
                R(j, j) = norm(v);
                Q(:, j) = v / R(j, j);
                %QRFactorization.plotQ(Q(:,j),mesh,bc,j)
                %Ap = Q(:,j).*R;
                Ap = Q(:,1:j)*R(1:j,1:j);
                x(1:j,1) = (Ap'*Ap)\(Ap'*RHS);
                res = Ap*x(1:j,1) - RHS;
                QRFactorization.plotRes(x,mesh,bc,j)
            end
        end

        function plotQ(Q,mesh,bc,iter)
            xFull = bc.reducedToFullVector(Q);
            s.fValues = reshape(xFull,2,[])';
            s.mesh = mesh;
            s.fValues(:,end+1) = 0;
            s.ndimf = 3;
            xF = P1Function(s);
            %xF.plot();
            xF.print(['Q',num2str(iter)],'Paraview')
            fclose('all');
        end

        function plotRes(res,mesh,bc,iter)
            xFull = bc.reducedToFullVector(res);
            s.fValues = reshape(xFull,2,[])';
            s.mesh = mesh;
            s.fValues(:,end+1) = 0;
            s.ndimf = 3;
            xF = P1Function(s);
            %xF.plot();
            xF.print(['Res',num2str(iter)],'Paraview')
            fclose('all');
        end

    end

end