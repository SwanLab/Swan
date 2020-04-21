function AnalyticalVsNumerical

sE = PonderatedOptimalSuperEllipseComputer();
sE.compute();

t = abs(sE.mxV - sE.myV) < 1e-6;
qV = sE.qMean(t);
m1 = sE.mxV(t);
m2 = sE.myV(t);

qA(:,1) = computeQanalytic(m1,m2);

f = figure();
hold on
h{1} = plot(m1,qA,'-+');
h{2} = plot(m1,qV,'-+');
xlabel('$m_1 = m_2$','Interpreter','Latex')
ylabel('$q$','Interpreter','Latex')

legA = '$\textrm{Analytical smoothing exponent} \, q_A$';
legN = '$\textrm{Numerical smoothing exponent} \, q_N$';

legObj = legend({legA,legN},'Interpreter','Latex','Location','Best');

outPutPath = '/home/alex/git-repos/MicroStructurePaper/';
outputName = [outPutPath,'qMaxM1M2AnalyticalNumerical'];
printer = plotPrinter(f,h);
printer.print(outputName);



end

function qA = computeQanalytic(mx,my)


for ipoint = 1:length(mx)
    
    s.m1 = mx(ipoint);
    s.m2 = my(ipoint);
    s.type = 'Optimal';
    qExp = SmoothingExponentComputer.create(s);
    qA(ipoint) = qExp.compute();
    
    
end

end