function AnalyticalVsNumerical

v = VademecumReader();
s.vademecum = v;

sE = PonderatedOptimalSuperEllipseComputer(s);
sE.compute();

%t = abs(v.mxV(:,1) - v.myV(:,1)) < 1e-6;

m1 = v.mxV(:,1);
m2 = v.myV(:,1);
q  = sE.qMean;

thet = 45;
theta = thet*pi/180;
m1L = @(x) sqrt(x^2*(1+tan(theta)^2));
m2L = @(y) sqrt(y^2*(1+1/(tan(theta)^2)));
tmin = min(m1L(max(m1)),m2L(max(m2)));
tmax = min(m1L(max(m1)),m2L(max(m2)));

t = linspace(0,1,100);
m1t = min(m1) + t/tmax*cos(theta);
m2t = min(m2) + t/tmax*sin(theta);

F = scatteredInterpolant(m1,m2,q);
qI = F(m1t,m2t);


qA(:,1) = computeQanalytic(m1t,m2t);

f = figure();
hold on
h{1} = plot(t,qA,'-+');
h{2} = plot(t,qI,'-+');
xlabel(['$\rho \quad (\xi = ',num2str(thet),') $'],'Interpreter','Latex')
ylabel('$q$','Interpreter','Latex')

legA = '$\textrm{Analytical smoothing exponent} \, q_A$';
legN = '$\textrm{Numerical smoothing exponent} \, q_N$';

legObj = legend({legA,legN},'Interpreter','Latex','Location','Best');

outPutPath = '/home/alex/git-repos/MicroStructurePaper/';
outputName = [outPutPath,'qMaxM1M2AnalyticalNumerical'];
printer = plotPrinter(f,h);
%printer.print(outputName);



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