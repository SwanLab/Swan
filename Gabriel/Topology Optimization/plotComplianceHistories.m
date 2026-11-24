% function plotComplianceHistories(result,fileName,normalizeBySolid)
function plotComplianceHistories(result,fileName,normalizeBySolid,volumeFraction)
if nargin < 2
    fileName = '';
end
if nargin < 3 || isempty(normalizeBySolid)
    normalizeBySolid = false;
end
if nargin < 4
    volumeFraction = [];
end
scale = 1;
yText = 'Compliance $\frac{1}{2}\int_\Omega \varepsilon : C : \varepsilon \, d\Omega$';
if normalizeBySolid
    scale = min(arrayfun(@(r) r.compliance(1),result));
    yText = 'Compliance / compliance of the solid';
end
figure('Color','w');
hold on
for k = 1:numel(result)
    J = result(k).compliance/scale;
    plot(0:numel(J)-1,J,'LineWidth',1.6, ...
        'DisplayName',sprintf('%s  (final %.4g, %d iterations)', ...
        result(k).label,J(end),numel(J)-1));
end
set(gca,'YScale','log');
grid on
xlabel('Iteration','Interpreter','latex');
ylabel(yText,'Interpreter','latex');
legend('Location','northeast','Interpreter','latex');
% title('Volume fraction 0.40','Interpreter','latex');
if ~isempty(volumeFraction)
    title(sprintf('Volume fraction %.2f',volumeFraction),'Interpreter','latex');
end
if ~isempty(fileName)
    print(gcf,fileName,'-dpng','-r300');
end
end