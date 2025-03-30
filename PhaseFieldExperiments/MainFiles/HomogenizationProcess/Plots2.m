clc,clear,close all
matType{1} = load('CircleMicroDamageArea.mat');
matType{2} = load('CircleMicroDamagePerimeter.mat');
matType{3} = load('SquareMicroDamageArea.mat');
matType{4} = load('SquareMicroDamagePerimeter.mat');
matType{5} = load('IsoMicroDamage.mat');
matType{6} = load('HorizontalCrackDamage.mat');
matType{6}.mat = matType{6}.mat*210;

% Change of variable
%matType{6}.phi  = matType{6}.holeParam{1};
%matType{6}.phi = matType{6}.holeParam{1}.^2;
matType{6}.phi = matType{6}.phi{1};

%% Polynomial fraction fitting
phi = matType{6}.holeParam{1};
C = squeeze(matType{6}.mat(3,3,:));

num = @(p) (  p(21).*phi.^10 + p(19).*phi.^9 + p(17).*phi.^8 + p(15).*phi.^7 + p(13).*phi.^6 ...
            + p(11).*phi.^5 + p(9).*phi.^4 + p(7).*phi.^3 + p(5).*phi.^2 ...
            + p(3).*phi.^1 + p(1));
den = @(p) (  p(22).*phi.^10 + p(20).*phi.^9 + p(18).*phi.^8 + p(16).*phi.^7 + p(14).*phi.^6 ...
            + p(12).*phi.^5 + p(10).*phi.^4 + p(8).*phi.^3 + p(6).*phi.^2 ...
            + p(4).*phi.^1 + p(2));
Cp = @(p) num(p)./den(p) + p(23);


objective = @(p) sum(sqrt(((Cp(p)'-C)./C).^2));
options = optimoptions(@fminunc,'StepTolerance',1e-10,'OptimalityTolerance',1e-10,...
                        'maxFunctionEvaluations',5e3,'MaxFunctionEvaluations',10000);

objResUnc = 100;
for i=1:1000
p0 = rand(1,23);
[popt,fval,exitflag,output,grad,hessian] = fminunc(objective,p0,options);
if fval<objResUnc
    objResUnc = fval;
    pResUnc = popt;
end
end

p = pResUnc;
num = @(x) (  p(21).*x.^10 + p(19).*x.^9 + p(17).*x.^8 + p(15).*x.^7 + p(13).*x.^6 ...
            + p(11).*x.^5 + p(9).*x.^4 + p(7).*x.^3 + p(5).*x.^2 ...
            + p(3).*x.^1 + p(1));
den = @(x) (  p(22).*x.^10 + p(20).*x.^9 + p(18).*x.^8 + p(16).*x.^7 + p(14).*x.^6 ...
            + p(12).*x.^5 + p(10).*x.^4 + p(8).*x.^3 + p(6).*x.^2 ...
            + p(4).*x.^1 + p(2));
Cp = @(x) num(x)./den(x) + p(23);

disp("Initial objective: " + num2str(objective(p0)));
disp("Final objective: " + num2str(objective(pResUnc)));
plot(phi,C,'ro')
hold on
fplot(Cp,[0 1],'-')
legend('measured','optimal')

%% Polynomial fraction fitting constrained
tiledlayout(3,3)
fun = cell(3,3);
for i=1:3
    for j=1:3
        phi = matType{6}.phi;
        C = squeeze(matType{6}.mat(i,j,:));

        num = @(p) (p(21).*phi.^10 + p(19).*phi.^9 + p(17).*phi.^8 + p(15).*phi.^7 + p(13).*phi.^6 ...
            + p(11).*phi.^5 + p(9).*phi.^4 + p(7).*phi.^3 + p(5).*phi.^2 ...
            + p(3).*phi.^1 + p(1));
        den = @(p) (p(22).*phi.^10 + p(20).*phi.^9 + p(18).*phi.^8 + p(16).*phi.^7 + p(14).*phi.^6 ...
            + p(12).*phi.^5 + p(10).*phi.^4 + p(8).*phi.^3 + p(6).*phi.^2 ...
            + p(4).*phi.^1 + p(2));
        Cp = @(p) num(p)./den(p);

        A = []; b = []; Aeq = []; beq = []; lb = []; ub = [];

        if i==1 && j==1
            nonlcon = @(p) nonLinearCon11(p,C);
        else
            nonlcon = @(p) nonLinearCon(p,C);
        end

        objective = @(p) sum(sqrt(((Cp(p)'-C)./C).^2));
        options = optimoptions(@fmincon,'StepTolerance',1e-10,'OptimalityTolerance',1e-10,...
            'ConstraintTolerance',1e-10,'MaxFunctionEvaluations',10000);

        objResCon = 100;
        for k=1:1000
            k
            p0 = rand(1,22);
            [popt,fval,exitflag,output,lambda,grad,hessian] = fmincon(objective,p0,A,b,Aeq,beq,lb,ub,nonlcon,options);
            if fval<objResCon
                objResCon = fval;
                pResCon = popt;
            end
        end

        p = pResCon;
        syms x
        num = (p(21).*x.^10 + p(19).*x.^9 + p(17).*x.^8 + p(15).*x.^7 + p(13).*x.^6 ...
            + p(11).*x.^5 + p(9).*x.^4 + p(7).*x.^3 + p(5).*x.^2 ...
            + p(3).*x.^1 + p(1));
        den = (p(22).*x.^10 + p(20).*x.^9 + p(18).*x.^8 + p(16).*x.^7 + p(14).*x.^6 ...
            + p(12).*x.^5 + p(10).*x.^4 + p(8).*x.^3 + p(6).*x.^2 ...
            + p(4).*x.^1 + p(2));
        Cp = num./den;
        fun{i,j} = Cp;

        disp("Initial objective: " + num2str(objective(p0)));
        disp("Final objective: " + num2str(objective(pResCon)));

        nexttile
        plot(phi,C,'X','LineWidth',1.5)
        hold on
        fplot(Cp,[0 1],'LineWidth',1.5)

        ylabel(char(8450)+"12 [GPa]");
        ylim([0,inf])
        xlabel("Damage "+char(632)+" [-]");
    end
end
lg =legend('measured','optimal');
lg.Layout.Tile = 'East';

dfun  = cell(3,3);
ddfun = cell(3,3);
for i=1:3
    for j=1:3
        dfun{i,j} = diff(fun{i,j});
        ddfun{i,j} = diff(dfun{i,j});
        if (i==1 && j==3) || (i==2 && j==3) || (i==3 && j==1) || (i==3 && j==2)
            degradationFun.fun{i,j} = @(x) x-x;
            degradationFun.dfun{i,j} = @(x) x-x;
            degradationFun.ddfun{i,j} = @(x) x-x;
        else
            degradationFun.fun{i,j} = matlabFunction(fun{i,j});
            degradationFun.dfun{i,j} = matlabFunction(dfun{i,j});
            degradationFun.ddfun{i,j} = matlabFunction(ddfun{i,j}); 
        end

    end
end
%% Include points
% %% Include final points
% matType{6}.mat(:,:,end+1) = [matType{6}.mat(1,1,end), 0 , 0;
%                                         0           , 0 , 0;
%                                         0           , 0 , 0];
% matType{6}.phi(end+1) = 1;
% matType{6}.mat(1,1,end) = matType{6}.mat(1,1,end-1);
%
% %% Include more final points %%
% lastSteps = 2;
% nPoints = 20;
% lastPointsMat = matType{6}.mat(:,:,(end-lastSteps):end);
% lastPointsPhi = matType{6}.phi((end-lastSteps):end);
% mat = zeros(3,3,lastSteps*(nPoints-1)+1);
% phi = zeros(1,lastSteps*(nPoints-1)+1);
% for k=1:lastSteps
%     if k==1 
%         extra = 1; 
%     else 
%         extra=0;
%     end
% 
%     for i=1:3
%         for j=1:3
%             start  = lastPointsMat(i,j,k);
%             finish = lastPointsMat(i,j,k+1);
%             res = linspace(start,finish,nPoints);
%             mat(i,j,((k-1)*nPoints+extra:k*nPoints-1)) = res(1:end-extra);
%         end
%     end
%     start  = lastPointsPhi(k);
%     finish = lastPointsPhi(k+1);
%     res = linspace(start,finish,nPoints);
%     phi(((k-1)*nPoints+extra:k*nPoints-1)) = res(1:end-extra);
% end
% mat(:,:,end) = lastPointsMat(:,:,end);
% phi(end) = lastPointsPhi(end);
% 
% matType{6}.mat = cat(3,matType{6}.mat(:,:,1:(end-lastSteps-1)),mat);
% matType{6}.phi = [matType{6}.phi(1:(end-lastSteps-1)),phi];

%% Fitting
 for j=1:length(matType)
     isMat6 = false;
     if j==6 isMat6 = true; end
     [funMat(:,:,j),dfunMat(:,:,j),ddfunMat(:,:,j)] = computeFittingHomogenization(matType{j},isMat6,4);
 end

 %% Run Plots
 Finalplots;


function [c,ceq] = nonLinearCon(p,y)
c = [];
ceq(1) = p(1)/p(2) - y(1);
ceq(2) = (p(1)+p(3)+p(5)+p(7)+p(9)+p(11)+p(13)+p(15)+p(17)+p(19)+p(21))/...
         (p(2)+p(4)+p(6)+p(8)+p(10)+p(12)+p(14)+p(16)+p(18)+p(20)+p(22));
end


function [c,ceq] = nonLinearCon11(p,y)
c = [];
ceq(1) = p(1)/p(2) - y(1);
ceq(2) = (p(1)+p(3)+p(5)+p(7)+p(9)+p(11)+p(13)+p(15)+p(17)+p(19)+p(21))/...
         (p(2)+p(4)+p(6)+p(8)+p(10)+p(12)+p(14)+p(16)+p(18)+p(20)+p(22)) - y(end);
end