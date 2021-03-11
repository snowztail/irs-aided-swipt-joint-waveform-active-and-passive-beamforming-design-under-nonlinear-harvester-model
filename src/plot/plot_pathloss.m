clear; clc; close all;

% number of distance samples
nSamples = 1e2;
% AP-user distance
directDistance = 12;
% vertical distance from the IRS to the AP-user path
verticalDistance = 2;
% projection of AP-IRS distance to the AP-user path
horizontalDistance = linspace(0, directDistance, nSamples);

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

%% * Auxiliary path loss
figure('name', 'Auxiliary path loss product vs AP-IRS horizontal distance', 'position', [0, 0, 500, 400]);
plotHandle = plot(horizontalDistance, pow2db(1 ./ auxiliaryPathloss));
grid on;
legend('$\Lambda_I\Lambda_R$', 'location', 'ne')
xlabel('AP-IRS horizontal distance [m]');
ylabel('Auxiliary path loss product [dB]');
box on;

savefig('../figures/path_loss.fig');
matlab2tikz('../../assets/path_loss.tex', 'extraaxisoptions', ['title style={font=\huge}, ' 'label style={font=\huge}, ' 'ticklabel style={font=\LARGE}, ' 'legend style={font=\LARGE}']);
