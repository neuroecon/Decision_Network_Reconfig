function community_detection_local_loo(dataDir,resultDir,roilist, mode,density,gamma,omega,subList,suffix)
% pipeline for multi-slice community detection, after high-pass filter correction, using the
% closest partition instead of group 'concensus'
% target at 118 rois involved in the 'task' and 'DMN' modules
% 4 slices: resting-state, semantic fluency, external-menu choice, internal-menu choice,
% density: network density
% gamma: spatial parameter
% omega: temporal parameter
% consthres: threshold of concensus partitioning
% Toolbox used: Brain Connectivity Toolbox, GenLovain, 
% Qianying Wu, 07/05/2025

condition={'Rest','SF','EMC','IMC'};
rsDir=[dataDir,condition{1},'/FCmat/'];
sfDir=[dataDir,condition{2},'/FCmat_',mode,'/'];
emcDir=[dataDir,condition{3},'/FCmat_',mode,'/'];
imcDir=[dataDir,condition{4},'/FCmat_',mode,'/'];
param=[num2str(density),',',num2str(gamma),',',num2str(omega)];

save2dir=[resultDir,'community_detection/LOO/',param,'_',mode,'no_',suffix,'/Subnetwork/'];
mkdir(save2dir);

%community detection
for sub=1:length(subList)
    % load FC matrix
    zmat{1}=load([rsDir,'zSub',int2str(subList(sub)),'.txt']); % Resting
    if strcmp(mode,'concat')
        zmat{2}=load([sfDir,'zSub',int2str(subList(sub)),'.txt']); % SF
        zmat{3}=load([emcDir,'zSub',int2str(subList(sub)),'.txt']); % EMC
        zmat{4}=load([imcDir,'zSub',int2str(subList(sub)),'.txt']); % IMC
    else % 'avg'
        zmat{2}=load([sfDir,'zSub',int2str(subList(sub)),'.mat']); % SF
        zmat{3}=load([emcDir,'zSub',int2str(subList(sub)),'.mat']); % EMC
        zmat{4}=load([imcDir,'zSub',int2str(subList(sub)),'.mat']); % IMC
    end        
        
    %reconstruct the network with 118 rois and threshold at given density
    for j=1:4
        zmat{j}=zmat{j}(roilist,:);
        zmat{j}=zmat{j}(:,roilist);
        zthres{j}=threshold_proportional(zmat{j},density);
    end
    
    % repeat 100 times
    for i=1:100
        [S,Q]=MultiCommunity(zthres,gamma,omega);
        Qallps(i)=Q;
        for j=1:4
            Sallps(:,i,j)=S(:,j); % each layer for 1 slice
        end
    end
    
    Qall(sub,:)=Qallps;
    
    % save Sallps
    savedir=[save2dir,'Sub',int2str(subList(sub)),'/'];
    if ~exist(savedir)
        mkdir(savedir);
    end
    save([savedir,'localS100.mat',],'Sallps');
    save([savedir,'zmatlocal.mat',],'zmat');
    
    Sall(:,sub,:)=minVInpartition(Sallps,4);
   
end
% save group level partition 
save([save2dir,'S27.mat'],'Sall');
save([save2dir,'Qall.mat'],'Qall');
    
%find a concensus across participants
ciu=minVInpartition(Sall,4);

save([save2dir,'partition_local.mat'],'ciu');