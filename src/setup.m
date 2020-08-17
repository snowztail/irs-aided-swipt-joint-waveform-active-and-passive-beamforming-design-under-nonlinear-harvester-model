addpath(genpath(pwd));
addpath('/rdsgpfs/general/user/yz16718/home/code/cvx');
cvx_setup /rdsgpfs/general/user/yz16718/home/code/cvx/cvx_license.dat
LASTN = maxNumCompThreads('automatic');
% rng(str2num(getenv('PBS_ARRAY_INDEX')));
