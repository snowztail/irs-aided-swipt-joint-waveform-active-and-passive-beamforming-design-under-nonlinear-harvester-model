addpath(genpath(pwd));
run('/rds/general/user/yz16718/home/code/cvx/cvx_setup.m');

nCpus = str2num(getenv('NCPUS'));
parpool(nCpus);
maxNumCompThreads(nCpus);
feature('numcores');
