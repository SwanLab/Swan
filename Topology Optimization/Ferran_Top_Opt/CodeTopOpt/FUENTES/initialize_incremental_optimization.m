function incropt = initialize_incremental_optimization(TYPE)

%% Perimeter
nsteps = 1;
incropt.alpha_eps = generate_incr_sequence(0,1,nsteps,'custom',8);  

%% C-C*
nsteps = 1;
if strcmp(TYPE,'MACRO')
    nsteps = 1;
end
incropt.alpha_Cstar = generate_incr_sequence(1/nsteps,1,nsteps,'custom',5);

%% Isotropy 2D
nsteps = 1;
if strcmp(TYPE,'MACRO')
    nsteps = 1;
end
incropt.alpha_isotropy2d = generate_incr_sequence(0,1,nsteps,'linear');

%% Volume
nsteps = 1;
incropt.alpha_vol = generate_incr_sequence(1/nsteps,1,nsteps,'linear');

%% Steps unification
varnames = fieldnames(incropt); % list of all incremental variables
fun = @(x0,x) [x0,length(x)];
nsteps = [];
for i = 1:length(varnames)
    eval(['nsteps = feval(fun,nsteps,incropt.',varnames{i},');']);
end
max_nsteps = max(nsteps);
for i = 1:length(varnames)
    eval(['incropt.',varnames{i},'(nsteps(i)+1:max_nsteps) = 1;']);
end
message = sprintf('Total number of steps = %g',max_nsteps);
try
    cprintf('err','%s\n',message);
catch ME
    warning('cprintf failed!');
    fprintf('%s\n',message);
end

%% Tolerances
incropt.alpha_PrimalTol = generate_incr_sequence(0,1,max_nsteps,'linear');
incropt.alpha_DualTol = generate_incr_sequence(0,1,max_nsteps,'linear');

end

function x = generate_incr_sequence (x1,x2,nsteps,type,factor)

switch type
    case 'linear'
        x = linspace(x1,x2,nsteps);
        
    case 'epsilon_sequence'
        frac = 2;
        kmax = ceil(log10(x1/x2)/log10(frac));
        x = epsilon0./frac.^(1:kmax);
        
    case 'logarithmic'
        x = logspace(log10(x1),log10(x2),nsteps);
        
    case 'custom'
        if nsteps < 2
            x = x2;
        else
            isteps = 0:nsteps-1;
            x = 1-(1-isteps/(nsteps-1)).^(factor);
            x = (x2-x1)*x + x1;
        end
        
    case 'free'
        x = zeros(1,nsteps);
        x(end) = 1;
        
    otherwise
        error('Incremental sequence type not detected.')
end

end