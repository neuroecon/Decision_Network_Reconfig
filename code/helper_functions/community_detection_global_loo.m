function community_detection_global_loo(dataDir,resultDir,mode, density,gamma,omega,subList,suffix)
% pipeline for global level multi-slice community detection 
% using the most similar partition as 'concensus'
% 4 slices: resting-state, semantic fluency, external-menu choice, internal-menu choice,
% dataDir: directory to read connectivity matrix
% resultDir: directory to save community detection output
% mode: FC matrix is concatenated ('concat') or averaged ('avg')
% density: network density
% gamma: spatial parameter
% omega: temporal parameter
% Toolbox used: Brain Connectivity Toolbox, GenLovain, 
% Qianying Wu, 07/05/2025

condition={'Rest','SF','EMC','IMC'};
rsDir=[dataDir,condition{1},'/FCmat/'];
sfDir=[dataDir,condition{2},'/FCmat_',mode,'/'];
emcDir=[dataDir,condition{3},'/FCmat_',mode,'/'];
imcDir=[dataDir,condition{4},'/FCmat_',mode,'/'];

param=[num2str(density),',',num2str(gamma),',',num2str(omega)];

save2dir=[resultDir,'community_detection/LOO/',param,'_',mode,'no_',suffix,'/'];
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
    
    %thresholding at given density
    for j=1:4
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
    savedir=[save2dir,'/Sub',int2str(subList(sub)),'/'];
    if ~exist(savedir)
        mkdir(savedir);
    end
    save([savedir,'S100.mat'],'Sallps');
    
    Sall(:,sub,:)=minVInpartition(Sallps,4);
   
end

% save group level partition
save([save2dir,'S27.mat'],'Sall');
save([save2dir,'Qall.mat'],'Qall');
    
%find a concensus across participants
ciu=minVInpartition(Sall,4);

save([save2dir,'partition.mat'],'ciu');

