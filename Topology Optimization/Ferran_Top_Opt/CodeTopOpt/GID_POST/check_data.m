function [ok] = check_data(selem)
ok=1;
if isfield(selem,'conec')
   else
   s=sprintf('table of connectivites is requered');
   disp(s); ok=0; 
   return
end
if isfield(selem,'matno')
   else
   s=sprintf('list of material number is requered');
   disp(s);ok=0; 
   return
end
if isfield(selem,'etype')
   else
   s=sprintf('element type is requered');
   disp(s);ok=0; 
   return
end
