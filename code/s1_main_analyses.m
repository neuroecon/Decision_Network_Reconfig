%% set the working directory as the codebook directory
wd = pwd;
addpath([wd,'/helper_functions']);
dataDir='../data/';
resultDir = '../results/';
condlist = {'Rest','SF','EMC','IMC'};
%% global level: network module detection 
% define parameters
density = 0.1;
gamma = 1;
omega = 0.1;
mode = 'concat';
community_detection_global(dataDir,resultDir,mode, density,gamma,omega)

%% global level: between-condition community similarity analysis

% 1. calculate Variation of Information (VIn) as a quantitative measure of
% network similairty
datadir = '../results/community_detection/0.1,1,0.1_concat';
load([datadir,'/S27.mat']);
subList=1:27;
for sub=1:length(subList)
    for i=1:3
        for j=i+1:4
            [VIn(i,j,sub),~] = partition_distance(Sall(:,sub,i), Sall(:,sub,j));
        end
    end
end
VIn_mean = mean(VIn, 3); % mean
VIn_std = std(VIn,[],3); % STD

% 2. perform Multi-Response Permutation Procedure to test the statistical
% significance of the differences
perm = 5000; % number of permutations
[rw0,rpwall,pval]=pipe_ANOS_perm(datadir,perm);

%% global level: network persistence analysis using the Lumped Markov Chain method
% need to first create ciu2 (manually), ciu2 is ciu after merging isolated
% nodes into an 'unknown' module (not an actual module), and assigning
% module numbers to the other modules
% 1=visual 2=FP Task 3=DMN 4=Sensory 5=Unknown
load('../results/community_detection/0.1,1,0.1_concat/partition.mat');
[mu_all, mu_mean, mu_SD]=pipe_Persistence_global(dataDir,ciu2,'concat');

%% module merge: network node flexibility
% all nodes for all subject (all partitions)
datadir = '../results/community_detection/0.1,1,0.1_concat';
subList=[702 705 708 711 718 719 720 722 725 728 729 730 733 735 ...
    736 737 740 744 745 748 750 754 755 758 759 760 761];
for sub=1:length(subList)
    load([datadir,'/Sub',int2str(subList(sub)),'/S100.mat']);
    for i=1:100 
        tempm=(reshape(Sallps(:,i,:),264,4))';
        F(:,i,sub) = flexibility(tempm,'cat');
    end 
end
Fmean1=reshape(mean(F,2),264,27);
Fgrandmean=mean(Fmean1,2);

% flexibility of modules (using the resting-state module partition)
load('../results/community_detection/0.1,1,0.1_concat/partition.mat');
F_visual = Fgrandmean(ciu2(:,1) == 1);
disp(['Flexibility of visual module = ',num2str(round(mean(F_visual),3)),'+-',num2str(round(std(F_visual),3))]);
F_FP = Fgrandmean(ciu2(:,1) == 2);
disp(['Flexibility of fronto-parietal module = ',num2str(round(mean(F_FP),3)),'+-',num2str(round(std(F_FP),3))]);
F_DMN = Fgrandmean(ciu2(:,1) == 3);
disp(['Flexibility of DMN module = ',num2str(round(mean(F_DMN),3)),'+-',num2str(round(std(F_DMN),3))]);
F_sensory = Fgrandmean(ciu2(:,1) == 4);
disp(['Flexibility of Sensory module = ',num2str(round(mean(F_sensory),3)),'+-',num2str(round(std(F_sensory),3))]);

%% module merge: nodal transition probability
% from RS FP to others in other conditions
load('../results/community_detection/0.1,1,0.1_concat/partition.mat');
jump =  [];
p_values = [];
for cond=2:4
    for mod=1:4
        jump(cond-1, mod) = sum(( ciu2(:,1) == 2) .* (ciu2(:,cond) == mod));
    end
end
jump_p = jump./38;

% perform Fisher's exact test comparing the transition probabilities
for cond=1:3
    for mod = [1,4]
        [~,p_values(cond,mod),~] = fishertest([jump(cond,3),38-jump(cond,3);jump(cond,mod),38-jump(cond,mod)]);
    end
end

%% module merge: compare the optimality of FP merge to DMN vs. merge to others
load('../results/community_detection/0.1,1,0.1_concat/partition.mat');

% combine task(2) and DMN(3)
ciu_task_dmn = ciu2;
for i =1:264
    if (sum(ciu2(i,1)==2) + sum(ciu2(i,1)==3)) > 0
        ciu_task_dmn(i,:) = [2,2,2,2];
    end
end
ciu_task_dmn(ciu_task_dmn == 5) = 3;
% persistence score for the merged module (row 2)
[mu_all, mu_mean, mu_SD]=pipe_Persistence(dataDir,ciu_task_dmn,'concat'); 

% combine task(2) and visual(1)
ciu_task_visual = ciu2;
for i =1:264
    if (sum(ciu2(i,1)==2) + sum(ciu2(i,1)==1)) > 0
        ciu_task_visual(i,:) = [2,2,2,2];
    end
end
ciu_task_visual(ciu_task_visual == 5) = 1;
% combine task(2) and sensory(4)
ciu_task_sensory = ciu2;
for i =1:264
    if (sum(ciu2(i,1)==2) + sum(ciu2(i,1)==4)) > 0
        ciu_task_sensory(i,:) = [2,2,2,2];
    end
end
ciu_task_sensory(ciu_task_sensory == 5) = 4;

% calculate modularity after each type of merge
Q_task_dmn = [];
Q_task_visual = [];
Q_task_sensory = [];
subList=[702 705 708 711 718 719 720 722 725 728 729 730 733 735 ...
736 737 740 744 745 748 750 754 755 758 759 760 761];
for sub=1:length(subList)
    zmat{1}=load([dataDir,'Rest/FCmat/zSub',int2str(subList(sub)),'.txt']); % Resting
    zmat{2}=load([dataDir,'SF/FCmat_concat/zSub',int2str(subList(sub)),'.txt']); % SF
    zmat{3}=load([dataDir,'EMC/FCmat_concat/zSub',int2str(subList(sub)),'.txt']); % EMC
    zmat{4}=load([dataDir,'IMC/FCmat_concat/zSub',int2str(subList(sub)),'.txt']); % IMC

    for j=1:4
        zthres{j}=threshold_proportional(zmat{j},0.1);
        [Q_task_dmn(sub,j)] = modularity_func(zthres{j},ciu_task_dmn(:,j));
        [Q_task_visual(sub,j)] = modularity_func(zthres{j},ciu_task_visual(:,j));
        [Q_task_sensory(sub,j)] = modularity_func(zthres{j},ciu_task_sensory(:,j));
    end
end
% print results
for cond = [2:4]
    disp([condlist{cond},'=================='])
    disp(['Q FP-DMN = ',num2str(round(mean(Q_task_dmn(:,cond)),3)),'+-',num2str(round(std(Q_task_dmn(:,cond)),3))]);
    disp(['Q FP-Visual = ',num2str(round(mean(Q_task_visual(:,cond)),3)),'+-',num2str(round(std(Q_task_visual(:,cond)),3))]);
    disp(['Q FP-DMN = ',num2str(round(mean(Q_task_sensory(:,cond)),3)),'+-',num2str(round(std(Q_task_sensory(:,cond)),3))]);
end

%% between module connectivity analysis: between visual and task-evoked module
% task-evoked = dmn + fp
load('../results/community_detection/0.1,1,0.1_concat/partition.mat');
[BetweenFCmean, BetweenFCsum, FCall]=between_mod_FC_1st(dataDir,ciu2,0.1, 1, [2,3]);
bar([1:4],mean(BetweenFCmean,1))                
hold on
errorbar([1:4],mean(BetweenFCmean,1),...
    std(BetweenFCmean,[],1)/sqrt(27),std(BetweenFCmean,[],1)/sqrt(27));   
xticklabels({'Rest','SF','EMC','IMC'});
hold off
ylabel('Between Visual and Task-evoked Mean Connectivity');
for i=[2,4]
    [~,p(i),~,stats{i}] = ttest(BetweenFCmean(:,3), BetweenFCmean(:,i));
end

% combine mean SBC condition
%meannonSBC = mean(BetweenFCmean(:,[2,4]),2);
%meanSBC = BetweenFCmean(:,3);
%boxplot([meanSBC, meannonSBC], 'Labels',{'SBC','non-SBC'});
%ylabel('Between Visual and Task-evoked Mean Connectivity');
%[tnonsbc, pnonsbc] = ttest(BetweenFCmean(:,3), meannonSBC);

%% local level: network module detection

load([resultDir,'community_detection/0.1,1,0.1_concat/Subnetwork/roilist.mat']);
community_detection_local(dataDir,resultDir,roilist, 'concat',0.1,1,0.1);

%% local level: between-condition community similarity analysis

% 1. calculate Variation of Information (VIn) as a quantitative measure of
% network similairty
datadir = '../results/community_detection/0.1,1,0.1_concat/Subnetwork/';
load([datadir,'S27.mat']);
subList=1:27;
for sub=1:length(subList)
    for i=1:3
        for j=i+1:4
            [VIn(i,j,sub),~] = partition_distance(Sall(:,sub,i), Sall(:,sub,j));
        end
    end
end
VIn_mean = mean(VIn, 3); % mean
VIn_std = std(VIn,[],3); % STD

% 2. perform Multi-Response Permutation Procedure to test the statistical
% significance of the differences
perm = 5000; % number of permutations
[rw0,rpwall,pval]=pipe_ANOS_perm(datadir,perm);


%% local level: network persistence analysis using the Lumped Markov Chain method
% need to first create ciu2 (manually), ciu2 is ciu after merging isolated
% nodes into an 'unknown' module (not an actual module), and assigning
% module numbers to the other modules
% 1,2,3,4 correspond to Mod 1,2,3,4
datadir = '../results/community_detection/0.1,1,0.1_concat/Subnetwork/';
load('../results/community_detection/0.1,1,0.1_concat/Subnetwork/partition_local.mat');
[Qdiag, meanQ, sdQ]=pipe_Persistence_local(datadir,ciu2);

%% between module connectivity analysis: between Mod 3 and Mod 4
% mod 3 is social + memory network, mod 4 is valuation network
% we hypothesize that the between-module connectivity is higher in IMC than
% SF
density = 0.1;
load('../results/community_detection/0.1,1,0.1_concat/Subnetwork/partition_local.mat');
[BetweenFCmean34, BetweenFCsum34, FCall34]= between_mod_FC_2nd(dataDir,ciu2,density, 3, 4);
BetweenFCmean34(isnan(BetweenFCmean34))=0;
bar([1:4],nanmean(BetweenFCmean34,1))                
hold on
errorbar([1:4],nanmean(BetweenFCmean34,1),...
    nanstd(BetweenFCmean34,[],1)/sqrt(27),nanstd(BetweenFCmean34,[],1)/sqrt(27));   
xticklabels({'Rest','SF','EMC','IMC'});
hold off
ylabel({'Between module connectivity:','Social/memory & Valuation network'})
for i=[1,2,3]
    [~,p(i),~,stats{i}] = ttest(BetweenFCmean34(:,4), BetweenFCmean34(:,i));
end


