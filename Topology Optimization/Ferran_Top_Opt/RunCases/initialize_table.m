function problem2run = initialize_table(casenum,casename,Vf,Ptarget,E_plus,E_minus,nu_plus,nu_minus,save_directory)

problem2run = table;

% Convert to vertical array
casenum = horz2vert (casenum);
casename = horz2vert (casename);
Vf = horz2vert (Vf);
Ptarget = horz2vert (Ptarget);
E_plus = horz2vert (E_plus);
E_minus = horz2vert (E_minus);
nu_plus = horz2vert (nu_plus);
nu_minus = horz2vert (nu_minus);
save_directory = horz2vert (save_directory);

% Get total number of cases
N = size(casenum,1);
N = max([N,size(casename,1)]);
N = max([N,size(Vf,1)]);
N = max([N,size(Ptarget,1)]);
N = max([N,size(E_plus,1)]);
N = max([N,size(E_minus,1)]);
N = max([N,size(nu_plus,1)]);
N = max([N,size(nu_minus,1)]);

% Complete data
if isnumeric(casenum)
    casenum = num2cell(casenum);
    fun = @(num) num2str(num);
    casenum = cellfun(fun,casenum,'UniformOutput',false);
end

complete = @(x) complete_data(x,N);

problem2run.casenum = complete(casenum);
problem2run.casename = complete (casename);
problem2run.Vf = complete (Vf);
problem2run.Ptarget = complete (Ptarget);
problem2run.E_plus = complete (E_plus);
problem2run.E_minus = complete (E_minus);
problem2run.nu_plus = complete (nu_plus);
problem2run.nu_minus = complete (nu_minus);

% Save directory
save_directory = complete (save_directory);
casenum = complete(casenum);
fun = @(str,num) sprintf('%s/',str); % sprintf('%s\\Case %s\\',str,num);
problem2run.save_directory = cellfun(fun,save_directory,casenum,'UniformOutput',false);

end

function x = horz2vert (x)

% x = x(:);
[row,col] = size(x);
if col > row
    x = x';
end

end

function x = complete_data(x,N)

M = size(x,1);
if M < N
    x(M:N,1) = x(M);
end

end