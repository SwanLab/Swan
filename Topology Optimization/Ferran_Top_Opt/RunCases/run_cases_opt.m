function run_cases_opt(method,kernel,algorithm,TYPE,problem_type,incropt,init_design_variable,perimeter_case,problem,top_opt_path,bpause,bnotification)

if ~isempty(top_opt_path)
    addpath(genpath(top_opt_path));
end

%% Generate input data
switch kernel
    case 'all'
        kernel={'PDE';'P1_kernel';'P0_kernel'};    
    case 'P0'
        kernel={'P0_kernel'};
    case 'P1'
        kernel={'P1_kernel'};
    otherwise
        kernel={kernel};
end

if ~iscell(algorithm)
    switch algorithm
        case 'all'
            algorithm={'Projected_gradient';'level_set';'MMA';'fmincon';'GCMMA';'IPOPT'};
        case 'PG'
            algorithm={'Projected_gradient'};
        case 'LS'
            algorithm={'level_set'};
        otherwise
            algorithm = {algorithm};
    end
end

switch method
    case 'all'
        method={'SIMP_ALL';'SIMP'};
    otherwise
        method = {method};
end

%% Call Opt Top
global check_wrong_ending
N=size(problem,1);
A=length(algorithm);
K=length(kernel);
M=length(method);
numtot=N*A*K*M;

fprintf('Starting %g cases\n',numtot);
% exitObj = onCleanup(@() exitfun(problem.save_directory{1},top_opt_path)); % execute exitfun on exit
for n = 1:N
    num=problem.casenum{n};
    save_directory = problem.save_directory{n};
    if ~exist(save_directory,'dir')
        mkdir(save_directory);
    end
    fid1=fopen([save_directory,'RunTimes.txt'],'at'); % 'wt'
    fid2=fopen([save_directory,'InfoCases.txt'],'at'); % 'wt'
    for m = 1:M
        for a = 1:A
            for k = 1:K
                % Input data              
                matprop.E_plus = problem.E_plus(n);
                matprop.E_minus = problem.E_minus(n);
                matprop.nu_plus = problem.nu_plus(n);
                matprop.nu_minus = problem.nu_minus(n);

                Algorithm = algorithm{a};
                kernel_case = kernel{k};
                method_case = method{m};

                Vfrac = problem.Vf(n);
                Ptar = problem.Ptarget(n);
                file_name = problem.casename{n};

                % Create directory for files
                file_info = sprintf('Case %s - %s %s %s %s - Vf %g Ptarget %g E+ %g E- %g nu+ %g nu- %g',num,file_name,method_case,Algorithm,kernel_case,Vfrac,Ptar,matprop.E_plus,matprop.E_minus,matprop.nu_plus,matprop.nu_minus);
                file = sprintf('Case %s - %s %s %s %s',num,file_name,method_case,Algorithm,kernel_case);                        
                file_write = sprintf('%s%s',save_directory,file);
                mkdir(file_write);
                figure_title = sprintf('Case %s - %s %s %s %s - Vf %g Ptarget %g',num,file_name,method_case,Algorithm,kernel_case,Vfrac,Ptar);
                
                fh  = figure('Name',figure_title);
                exitObj = onCleanup(@() exitfun(save_directory,top_opt_path,fh)); % execute exitfun on exit
                mp = get(0, 'MonitorPositions');
                select_screen = 1;
                if size(mp,1) < select_screen
                    select_screen = size(mp,1);
                end
                width = mp(1,3);
                height = mp(1,4);
                size_screen_offset = round([0.005*width,0.07*height,-0.25*width,-0.15*height],0);
                set(fh,'Position',mp(select_screen,:) + size_screen_offset);
                fprintf(fid2,'%s\n',file_info);

                % Call Main Code
                tic
                check_wrong_ending = true;
                MAIN_OPT3 (file_write,file_name,Vfrac,Ptar,method_case,Algorithm,kernel_case,matprop,TYPE,problem_type,incropt,init_design_variable,perimeter_case,fh,select_screen);
                time = toc;
                fprintf(fid1,'Case %s - %g sec\n',num,time);

                % Save figure with results
                savefig(sprintf('%s.fig',file_write));  
                check_wrong_ending = false;
                
                % Show case
                fprintf('%s - %s - %s - %s\n',file_name,method_case,Algorithm,kernel_case);
                if bpause
                    pause;
                end
                if bnotification
                    load handel
                    sound(y,Fs)
                end
                close(fh);
                if N > 1
                    close all;
                end
            end
        end
    end
    fclose(fid1);
    fclose(fid2);
end


end

function exitfun(path,top_opt_path,fh)

global check_wrong_ending
if ~check_wrong_ending
    return;
end
    
% Executes when exiting function (due to normal ending or Ctrl+C)
status = fclose('all');
if status == 0
    fprintf('Files closed successfully!\n\n');
else
    fprintf('Unable to close files.\n\n');
end

% Check if abnormal ending
if isgraphics(fh)  
    % Get figure name
    fid = fopen([path,'InfoCases.txt'],'r');
    data = textscan(fid,'%s','Delimiter','\n');
    fname = data{:}{end};
    fclose(fid);
    pattern = ' - Vf [^]*';
    fname = regexprep(fname,pattern,'');
    savefig(fh,[path,fname,'.fig']);
    winopen(path);
    close(fh);
    
    % Info messages
    fprintf('Forced exit of the program, figure was saved and files closed.\n');
    fprintf('Remember:\n');
    fprintf('\t-The problem was encountered during execution of case: %s\n',fname);
    fprintf('\t-Add data in RunTimes.txt manually for interrupted case\n');
    fprintf('\t-Delete data in InfoCases.txt if necessary\n');
    fprintf('\n\n');
    
    % Complete RunTimes
    fid=fopen([path,'RunTimes.txt'],'at');
    time = toc;
    strpattern = 'Case %s - %*s %*s %*s %*s';
    num = textscan(fname,strpattern);
    fprintf(fid,'Case %s - %g sec\n',num{:}{:},time);
    fclose(fid);
    
    % Open files to edit in case of failed run
    %edit([path,'InfoCases.txt']);
    %edit([path,'RunTimes.txt']);
end

% Return to default path
if ~isempty(top_opt_path)
    rmpath(genpath(top_opt_path));
end

% clear global;

end

