addpath(genpath(pwd));

set(groot, 'defaulttextinterpreter', 'latex');
set(groot, 'defaultaxesticklabelinterpreter', 'latex');
set(groot, 'defaultaxesfontname', 'latex');
set(groot, 'defaultlegendinterpreter', 'latex');
set(groot, 'defaultlinelinewidth', 2);
set(groot, 'defaultstemlinewidth', 2);

% run('/rds/general/user/yz16718/home/cvx/cvx_setup.m');
% % run('/rds/general/user/yx2718/home/yang/cvx/cvx_setup.m');
% 
% iBatch = str2double(getenv('PBS_ARRAY_INDEX'));
% disp(iBatch);
% rng(iBatch);
