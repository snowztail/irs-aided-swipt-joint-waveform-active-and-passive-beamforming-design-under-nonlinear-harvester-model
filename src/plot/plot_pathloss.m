clear; clc; close all;

% number of distance samples
nSamples = 1e2;
% AP-user distance
directDistance = 15;
% vertical distance from the IRS to the AP-user path
verticalDistance = 2;
% projection of AP-IRS distance to the AP-user path
horizontalDistance = linspace(0, directDistance, nSamples);

% AP-IRS pathloss
incidentPathloss = zeros(nSamples, 1);
% production of AP-IRS and IRS-user pathloss
extraPathloss = zeros(nSamples, 1);
for iSample = 1 : nSamples
	% AP-IRS and IRS-user distance
	[incidentDistance, reflectiveDistance] = coordinate(directDistance, verticalDistance, horizontalDistance(iSample));
	incidentPathloss(iSample) = path_loss(incidentDistance);
	extraPathloss(iSample) = path_loss(incidentDistance) * path_loss(reflectiveDistance);
end

%% * Waveform amplitude
figure('name', 'Path loss vs AP-IRS horizontal distance');
pathlossPlot = tiledlayout(2, 1, 'tilespacing', 'compact');

% * Incident pathloss
nexttile;
semilogy(horizontalDistance, incidentPathloss);
grid on;
legend('IEEE TGn D')
xlabel('Distance [m]');
ylabel('Path loss');

% * Extra pathloss
nexttile;
semilogy(horizontalDistance, extraPathloss);
grid on;
legend('$\Lambda_I\Lambda_R$', 'location', 'se')
xlabel('AP-IRS horizontal distance [m]');
ylabel('Path loss product');

savefig('../figures/path_loss.fig');
matlab2tikz('../../assets/path_loss.tex');
