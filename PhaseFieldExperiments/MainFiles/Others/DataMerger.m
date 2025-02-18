s.reaction = [];
s.displacement.value = [];
s.damage.maxValue = [];
s.damage.field = [];
s.energy.extWork = [];
s.energy.intE = [];
s.energy.localDis = [];
s.energy.regDis = [];
s.iter.u = [];
s.iter.phi = [];
s.iter.stag = [];
s.cost = [];

for i=1:length(data)
    s.reaction = [s.reaction data{i}.reaction];
    s.displacement.value = [s.displacement.value data{i}.displacement.value];
    s.damage.maxValue = [s.damage.maxValue data{i}.damage.maxValue];
    s.damage.field = data{i}.damage.field;
    s.energy.extWork = [s.energy.extWork data{i}.energy.extWork];
    s.energy.intE = [s.energy.intE data{i}.energy.intE];
    s.energy.localDis = [s.energy.localDis data{i}.energy.localDis];
    s.energy.regDis = [s.energy.regDis data{i}.energy.regDis];
    s.iter.u = [s.iter.u data{i}.iter.u];
    s.iter.phi = [s.iter.phi data{i}.iter.phi];
    s.iter.stag = [s.iter.stag data{i}.iter.stag];
    s.cost = [s.cost data{i}.cost];
end

PhaseFieldPlotter(s)