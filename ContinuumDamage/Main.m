close all
tCD= TestingContinuumDamage();
data = tCD.data;

figure()
plot(data.displacement.value,data.damage.maxValue)

figure()
plot(data.displacement.value,data.q.maxValue)

figure()
plot(data.displacement.value,data.reaction)