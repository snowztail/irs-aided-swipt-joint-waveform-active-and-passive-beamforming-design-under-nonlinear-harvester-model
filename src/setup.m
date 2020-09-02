addpath(genpath(pwd));
run('/rds/general/user/yz16718/home/code/cvx/cvx_setup.m');

iBatch = str2double(getenv('PBS_ARRAY_INDEX'));
disp(iBatch);
rng(iBatch);
