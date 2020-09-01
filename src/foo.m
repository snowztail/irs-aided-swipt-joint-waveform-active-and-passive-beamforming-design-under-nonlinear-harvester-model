addpath(genpath(pwd));
run('/rds/general/user/yz16718/home/code/cvx/cvx_setup.m');

pbsIndex = str2double(getenv('PBS_ARRAY_INDEX'));
disp(pbsIndex);
rng(pbsIndex);

load('data/foo.mat');
reSet(:, pbsIndex) = pbsIndex;
save('data/foo.mat', 'reSet', '-append');
