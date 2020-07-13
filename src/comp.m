clear; clc; setup; config_comp; load('ref.mat'); load('gp_sdr.mat');

% ni_gp;
ni_sdr;

plot(niGpSample(1, :) / nSubbands, 1e6 * niGpSample(2, :), 'b-');
hold on;
plot(niSdrSample(1, :) / nSubbands, 1e6 * niSdrSample(2, :), 'r-');
hold on;
plot(reRef(1, :), 1e6 * reRef(2, :), 'k--');
