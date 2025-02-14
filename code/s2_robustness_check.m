% robustness check: changing density and FC construction methods
%% global level: network module detection
dataDir='../data/';
resultDir = '../results/';
gamma = 1;
omega = 0.1;
% averaged FC
for density = [0.05,0.1,0.2]
    community_detection(dataDir,resultDir,'avg',density,gamma,omega);
end
% concatenated FC
for density = [0.05,0.2]
    community_detection(dataDir,resultDir,'concat',density,gamma,omega);
end

%% global level: between-condition community similarity analysis
datadir = '../results/community_detection/';
perm = 5000; % number of permutations
gamma = 1;
omega = 0.1;
% averaged FC
for density = [0.05,0.1,0.2]
    params = [num2str(density),',',num2str(gamma),',',num2str(omega),'_avg'];
    [rw0,rpwall,pval]=pipe_ANOS_perm([datadir,params],perm);
end
% concatenated FC
for density = [0.05,0.2]
    params = [num2str(density),',',num2str(gamma),',',num2str(omega),'_concat'];
    [rw0,rpwall,pval]=pipe_ANOS_perm([datadir,params],perm);
end
