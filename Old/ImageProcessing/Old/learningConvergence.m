function learningConvergence

k = 10;
iter = 1:200;

lblRate  = LowerBoundLipshitzFirstOrder(iter);
lbscRate = lowerBoundStronglyConvexFirstOrder(iter,k);
sRate = subgradientConvergence(iter);
nRate = stronglyConvexNesterovConvergence(iter,k);
qRate = quadraticConvergence(iter,1/k);
lgRate = LipschitzGradientConvergence(iter);
scRate = stronglyConvexGradientConvergence(iter,k);

colors =    {
            [0.9290, 0.6940, 0.1250];	          	
            [0, 0.4470, 0.7410];
          	[0.8500, 0.3250, 0.0980];
            [0, 0.4470, 0.7410];
            [0.8500, 0.3250, 0.0980];
          	[0.4940, 0.1840, 0.5560];	          	
          	[0.4660, 0.6740, 0.1880];	          	
          	[0.3010, 0.7450, 0.9330]};
style = {'-';'--';'--';'-';'-';'-';'-';};       

rates = [sRate;lblRate;lbscRate;lgRate;scRate;nRate;qRate];

figure(1)
subplot(1,2,1);
hold on
for i = 1:size(rates,1)
    plot(iter,rates(i,:),'Color',colors{i},'LineStyle',style{i});
end
legend('Subgradient','LowerBound Lipshitz','LowerBound Strongly Convex', 'Lipschitz Diff Gradient','StronglyConvex Diff Gradient','Strongly Convex Diff Nesterov','Quadratic Diff')

figure(1)
subplot(1,2,2);
hold on
for i = 1:size(rates,1)
h = plot(iter,rates(i,:),'Color',colors{i},'LineStyle',style{i});
set(gca,'yscale','log');
end


legend('Subgradient','LowerBound Lipshitz','LowerBound Strongly Convex', 'Lipschitz Diff Gradient','StronglyConvex Diff Gradient','Strongly Convex Diff Nesterov','Quadratic Diff')


end

function rate = subgradientConvergence(x)

rate = 1./sqrt(x);

end

function rate = LowerBoundLipshitzFirstOrder(x)

rate = 1./(x + 1).^2;

end



function rate = LipschitzGradientConvergence(x)

rate = 1./(x + 1);

end

function rate = stronglyConvexGradientConvergence(x,k)

rate = (((k) - 1)/((k)+1)).^(x);
end


function rate = quadraticConvergence(x,k)

rate = k.^(2.^x);

end


function rate = lowerBoundStronglyConvexFirstOrder(x,k)

rate = ((sqrt(k) - 1)/(sqrt(k) + 1)).^(2*x);

end


function rate = stronglyConvexNesterovConvergence(x,k)

rate = ((sqrt(k) - 1)/(sqrt(k)+1)).^(x);

end