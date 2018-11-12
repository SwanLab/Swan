function readVademecum
Path = '/media/Alex/Vademecum';
nlevel = 7;
pint = 0; % 0 == no pintar, 1 == pintar
[~,~,~,numsna] = vademecum_spehere_mesh(nlevel,pint);

nsnapshots = numsna(end);

i_not_done = 1;
i_not_finished = 1;
i_not_converged = 1;
i_done = 1;
i_conv = 1;
case_not_done = [];
case_not_converged = [];
case_not_finished = [];
for isnap = 1:nsnapshots
    
    
    fid = fopen([Path,'/',num2str(isnap),'/executed_cases.txt'],'r');     %# Open the file as a binary
    
    if fid == -1
        case_not_done(i_not_done) = isnap;
        i_not_done = i_not_done + 1;
    else
        lastline = call_last_line(fid);
        fclose(fid);  %# Close the file
        C = textscan(lastline, '%s%s%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f');
        k = 1;
        for i = 2:2:28
            Datos{k} = [C{i}];
            k = k + 1;
        end
        %NewData = reshape(Datos,[],5);
        tipo = [Datos{1}];
        Phi = Datos{2};
        Theta = Datos{3};
        Vol = Datos{4};
        Lagrange0 = Datos{5};
        Iter = Datos{6};
        LagrangeInf = Datos{7};
        Cost0 = Datos{8};
        Cost_Inf = Datos{9};
        h_inf = Datos{10};
        theta_0 = Datos{11};
        theta_inf = Datos{12};
        Energia = Datos{13};
        Mesh = Datos{14};
        
        if isequal(num2str(Mesh),'NaN') || length(Mesh) == 0
            case_not_finished(i_not_finished)  = isnap;
            i_not_finished = i_not_finished + 1;
            %unix(['rm -r ',Path,'/',num2str(isnap)])
        else
            isnap
            thet(i_done) = theta_inf;
            ener(i_done) = Energia;
            phi_snap(i_done) = Phi;
            theta_snap(i_done) = Theta;
            iteracio(i_done) = isnap;
            i_done = i_done + 1;
            if theta_inf >= 4
                case_not_converged(i_not_converged) = isnap;
                i_not_converged = i_not_converged + 1;
               %unix(['rm -r ',Path,'/',num2str(isnap)])
            else
               converg(i_conv) = isnap;
               i_conv = i_conv + 1;
            end
            
        end
        
        
        
    end
    
end
to_do = [length(case_not_finished) length(case_not_done) length(case_not_converged)];
[num2str((nsnapshots-sum(to_do))/nsnapshots*100),'%']

end


function lastLine = call_last_line(fid)

lastLine = '';                   %# Initialize to empty
offset = 1;                      %# Offset from the end of file
fseek(fid,-offset,'eof');        %# Seek to the file end, minus the offset
newChar = fread(fid,1,'*char');  %# Read one character
while (~strcmp(newChar,char(10))) || (offset == 1)
    lastLine = [newChar lastLine];   %# Add the character to a string
    offset = offset+1;
    fseek(fid,-offset,'eof');        %# Seek to the file end, minus the offset
    newChar = fread(fid,1,'*char');  %# Read one character
end

end