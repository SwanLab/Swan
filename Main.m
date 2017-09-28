function obj = Main



%% Steps
% 1 - Run 'Main.m'
% 2 - Create object --> obj = Physical_Problem();
% 3 - Preprocess    --> obj.preProcess(filename_in);
% 4 - Compute       --> obj.computeVariables();
% 5 - Postprocess   --> obj.postProcess(filename_out);


%% Main.m

name_in  = 'CantileverToy_Tetrahedra';
name_out = 'results-toyExample'; 


folder_in  = 'Input';
folder_out = 'Output';

dir_in  = fullfile(pwd,folder_in);
dir_out = fullfile(pwd,folder_out);

filename_in  = fullfile(dir_in,name_in);
filename_out = fullfile(dir_out,name_out); 



obj = Physical_Problem();
obj.preProcess(filename_in);
obj.computeVariables();
obj.postProcess(filename_out);

end

