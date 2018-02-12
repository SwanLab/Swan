
function [conectivities,coordgid] = call_conectivities_1_over_8_sphere
eval('VademecumMesh')
conectivities = conectivities(:,2:end);
coordgid = coordinates;
end