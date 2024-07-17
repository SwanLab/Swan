clear;
close all;

fig = open("APReview\DensityPure.fig");

ch     = fig.Children;
trials = ch(4).Children.YData;
close all;

trialsUpdated = trials + 1;
nPDE = sum(trialsUpdated)+1;