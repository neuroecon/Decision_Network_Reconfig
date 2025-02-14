% between-module connectivity strength - global level
%1=visual 3=dmn 2=task 4=sensory 5=uncertain
function [BetweenFCmean, BetweenFCsum, FCall]= between_mod_FC_1st(dataDir,ciu2,density, mod1, mod2)

subList=[702 705 708 711 718 719 720 722 725 728 729 730 733 735 ...
736 737 740 744 745 748 750 754 755 758 759 760 761];
for sub=1:length(subList)
    zmat{1}=load([dataDir,'Rest/FCmat/zSub',int2str(subList(sub)),'.txt']); % Resting
    zmat{2}=load([dataDir,'SF/FCmat_concat/zSub',int2str(subList(sub)),'.txt']); % SF
    zmat{3}=load([dataDir,'EMC/FCmat_concat/zSub',int2str(subList(sub)),'.txt']); % EMC
    zmat{4}=load([dataDir,'IMC/FCmat_concat/zSub',int2str(subList(sub)),'.txt']); % IMC

    for j=1:4
        zthres{j}=threshold_proportional(zmat{j},density);
        [BetweenFCmean(sub,j),BetweenFCsum(sub,j), FCall{sub,j}]=BMstrength(zthres{j},ciu2(:,j),mod1,mod2); % between-module strength
    end
end



