
function [NMI] =run_hiC_network(ii,alpha_NR)


addpath(genpath(pwd))


HiC_raw_resolution1K = load(['mat_', num2str(ii),'_1000.txt']); % raw HiC contact network with resolution 1K
label_resolution1K = load(['class_',num2str(ii),'_1000.txt']); % labels for each node in HiC network with resolution 1K

HiC_NR=NR(HiC_raw_resolution1K,alpha_NR);

HiC_NR=Network_Enhancement(HiC_NR);

com_NR = louvain(HiC_NR);

NMI = nmi(com_NR,label_resolution1K);






