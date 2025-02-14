function [Qdiag, meanQ, sdQ]=pipe_Persistence_global(rootpath,ciu2,mode)
% Test the significance of individual community
% proposed by Carlo Piccardi(2011)'Finding and Testing Network Communities
% by Lumped Markov Chains'
% use alpha values to describe the significance, alpha threshold=0.5
% Qianying Wu, 02/01/2020

subList=[702 705 708 711 718 719 720 722 725 728 729 730 733 735 ...
736 737 740 744 745 748 750 754 755 758 759 760 761];

% 1=visual 2=FP Task 3=DMN 4=Sensory 5=Unknown
dataDir{1}=[rootpath,'Rest/FCmat/'];
dataDir{2}=[rootpath,'SF/FCmat_',mode,'/'];
dataDir{3}=[rootpath,'EMC/FCmat_',mode,'/'];
dataDir{4}=[rootpath,'IMC/FCmat_',mode,'/'];

for sub=1:length(subList)
    for cond=1:4
        clear z0 C
        C=ciu2(:,cond);
        if cond == 1
            z0=load([dataDir{cond},'zSub',int2str(subList(sub)),'.txt']); % load z-FC matrix
        else
            if strcmp(mode,'concat')
                z0 = load([dataDir{cond},'zSub',int2str(subList(sub)),'.txt']);
            else    
                z0 = load([dataDir{cond},'zSub',int2str(subList(sub)),'.mat']).FC_mat_avg; % load z-FC matrix
            end
        end
        z=threshold_proportional(z0,0.1);
        cdelete=find(all(z==0,2)==1);
        
        z(all(z==0,2),:)=[];
        z(:,all(z==0,1))=[]; % delete zero columns and rows
        
        C(cdelete,:)=[]; % delete zero columns from the partitions
        
        [H{cond,sub},Q{cond,sub},persistence{cond,sub}]=Sig_LMC(C,z);
        for col=1:size(Q{cond,sub},1)
            Qdiag{cond}(sub,col)=Q{cond,sub}(col,col);
        end
    end
end

% calculate mean and sd
meanQ=[];
sdQ=[];
for qcond=1:size(Qdiag,2)
    tempQ=Qdiag{qcond};
    for qpart=1:size(tempQ,2)
        tempartQ=tempQ(:,qpart);
        tempartQ=tempartQ(tempartQ<=1);
        tempartQ=tempartQ(tempartQ>=0);
        meanQ(qpart,qcond)=mean(tempartQ(tempartQ<=1));
        sdQ(qpart,qcond)=std(tempartQ(tempartQ<=1),0,1);
    end
end