function postprocess_macros(Path_macros,macros_evol,file_name,case_post,mirroring)
    if mirroring 
        tcl_file = ['createGidFiguresMirroring.tcl'];            
    else
        tcl_file = ['createGidFigures.tcl'];            
    end

    for imacro = 1:size(macros_evol,1)
        file_tcl_name = [Path_macros,'tcl_Video_file_macro_',num2str(imacro),'.tcl'];
        write_tcl_file(file_tcl_name,Path_macros,[Path_macros,file_name],size(macros_evol,2),case_post,tcl_file,macros_evol)
%       unix(['/opt/GiDx64/12.1.2d/gid_offscreen -t "source ',file_tcl_name,'"'])
        unix(['cp /home/aferrer/Documents/Doctorat/Tesi/FEM_DT/',tcl_file,' ',Path_macros])
        unix(['/opt/GiDx64/12.1.2d/gid_offscreen -offscreen -t "source ',file_tcl_name,'"']);
        crop_images_macro(macros_evol,[Path_macros,file_name]);
        unix(['rm ',file_tcl_name]);
    end
end

function write_tcl_file(file_tcl_name,Path,file,nstep,case_post,tcl_file,macros_evol)

file_macro = [];
for iter = 1:nstep
   file_macro = [file_macro,'"',file,'_',num2str(macros_evol(iter)),'.flavia.res" '];
end

fid = fopen(file_tcl_name,'w+');
fprintf(fid,'GiD_Process PostProcess \n');
fprintf(fid,['set Path_video_output "',file,'"\n']);
fprintf(fid,['set Path_macros_input [list ',file_macro,']\n']);

fprintf(fid,['set Name_figure ','macro','\n']);

switch case_post
    case 'video'    
        fprintf(fid,['source "',Path,'createGidVideo.tcl"\n']);
    case 'images'
        fprintf(fid,['source "',Path,tcl_file,'"\n']);            
end

deformation = 0;
switch file
    case [Path,'Perfil_Alar']
     deformation = 0.0005;   
    otherwise
     deformation = 0.023;   
end

fprintf(fid,['VideoPrueba ',num2str(deformation) ,' $Path_video_output $Path_macros_input $Name_figure\n']);
fprintf(fid,['GiD_Process Mescape Quit']);
fclose(fid);

end

function crop_images_macro(macros_evol,Path)
for istep = 1:size(macros_evol,2)
    name_file = [Path,'-',num2str(istep),'_macro','.png'];
    %unix(['convert -crop 1500x400+0+100 -gravity Center ',name_file,' ',name_file]);
    %unix(['convert -crop 1500x200+0+100 -gravity Center ',name_file,' ',name_file]);
    %unix(['convert -crop 700x200+0+0 -gravity Center ',name_file,' ',name_file]);
    unix(['convert -crop 800x300+0+0 -gravity Center ',name_file,' ',name_file]);
end

end