clear; clc; close all;

% number of distance samples
nSamples = 1e2;
% AP-user distance
directDistance = 12;
% vertical distance from the IRS to the AP-user path
verticalDistance = 2;
% projection of AP-IRS distance to the AP-user path
horizontalDistance = linspace(0, directDistance, nSamples);
% number of cases to plot (path loss, path loss product)
nCases = 2;

% AP-IRS pathloss
incidentPathloss = zeros(nSamples, 1);
% product AP-IRS-User pathloss model
auxiliaryPathloss = zeros(nSamples, 1);
for iSample = 1 : nSamples
	% AP-IRS and IRS-user distance
	[incidentDistance, reflectiveDistance] = coordinate(directDistance, verticalDistance, horizontalDistance(iSample));
	incidentPathloss(iSample) = path_loss(incidentDistance);
	auxiliaryPathloss(iSample) = path_loss(incidentDistance) * path_loss(reflectiveDistance);
end

%% * Waveform amplitude
figure('name', 'Path loss vs AP-IRS horizontal distance');
pathlossPlot = tiledlayout(nCases, 1, 'tilespacing', 'compact');
plotHandle = gobjects(1, nCases);

% * Incident pathloss
nexttile;
plotHandle(1) = plot(horizontalDistance, pow2db(1 ./ incidentPathloss));
grid on;
legend('$\Lambda_D$', 'location', 'se')
xlabel('AP-user distance [m]');
ylabel('Path loss [dB]');
box on;

% * Extra pathloss
nexttile;
plotHandle(2) = plot(horizontalDistance, pow2db(1 ./ auxiliaryPathloss));
grid on;
legend('$\Lambda_I\Lambda_R$', 'location', 'ne')
xlabel('AP-IRS horizontal distance [m]');
ylabel('Path loss product [dB]');
box on;

savefig('../figures/path_loss.fig');
matlab2tikz('../../assets/path_loss.tex');
