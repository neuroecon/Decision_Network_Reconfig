% pipeline of anosperm
function [rw0,rpwall,pval]=pipe_ANOS_perm(datadir,perm)

% rw0: actual stat
% rpall: all the stats after permutation
% pval: p-value

[rw12,rpw12,pvalue12]=ANOS_perm(datadir,1,2,perm);
[rw13,rpw13,pvalue13]=ANOS_perm(datadir,1,3,perm);
[rw14,rpw14,pvalue14]=ANOS_perm(datadir,1,4,perm);
[rw23,rpw23,pvalue23]=ANOS_perm(datadir,2,3,perm);
[rw24,rpw24,pvalue24]=ANOS_perm(datadir,2,4,perm);
[rw34,rpw34,pvalue34]=ANOS_perm(datadir,3,4,perm);
n=27;
subplot(2,3,1)
h1=histogram(rpw12,floor(perm/20),'EdgeColor',[0.3,0.75,0.93],'FaceColor',[0.3,0.75,0.93]);
hold on
plot(rw12*ones(n+1,1),linspace(0,max(h1.BinCounts)+3,n+1),'DisplayName','delta','LineWidth',1,...
    'Color',[0.635294139385223 0.0784313753247261 0.184313729405403]);
title('Rest vs SF');
hold on

subplot(2,3,2)
h2=histogram(rpw13,floor(perm/20),'EdgeColor',[0.3,0.75,0.93],'FaceColor',[0.3,0.75,0.93]);
hold on
plot(rw13*ones(n+1,1),linspace(0,max(h2.BinCounts)+3,n+1),'DisplayName','delta','LineWidth',1,...
    'Color',[0.635294139385223 0.0784313753247261 0.184313729405403]);
title('Rest vs EMC');
hold on

subplot(2,3,3)
h3=histogram(rpw14,floor(perm/20),'EdgeColor',[0.3,0.75,0.93],'FaceColor',[0.3,0.75,0.93]);
hold on
plot(rw14*ones(n+1,1),linspace(0,max(h3.BinCounts)+3,n+1),'DisplayName','delta','LineWidth',1,...
    'Color',[0.635294139385223 0.0784313753247261 0.184313729405403]);
title('Rest vs IMC');
hold on

subplot(2,3,4)
h4=histogram(rpw23,floor(perm/20),'EdgeColor',[0.3,0.75,0.93],'FaceColor',[0.3,0.75,0.93]);
hold on
plot(rw23*ones(n+1,1),linspace(0,max(h4.BinCounts)+3,n+1),'DisplayName','delta','LineWidth',1,...
    'Color',[0.635294139385223 0.0784313753247261 0.184313729405403]);
title('SF vs EMC');
hold on

subplot(2,3,5)
h5=histogram(rpw24,floor(perm/20),'EdgeColor',[0.3,0.75,0.93],'FaceColor',[0.3,0.75,0.93]);
hold on
plot(rw24*ones(n+1,1),linspace(0,max(h5.BinCounts)+3,n+1),'DisplayName','delta','LineWidth',1,...
    'Color',[0.635294139385223 0.0784313753247261 0.184313729405403]);
title('SF vs IMC');
hold on

subplot(2,3,6)
h6=histogram(rpw34,floor(perm/20),'EdgeColor',[0.3,0.75,0.93],'FaceColor',[0.3,0.75,0.93]);
hold on
plot(rw34*ones(n+1,1),linspace(0,max(h6.BinCounts)+3,n+1),'DisplayName','delta','LineWidth',1,...
    'Color',[0.635294139385223 0.0784313753247261 0.184313729405403]);
title('EMC vs IMC');

suptitle(datadir);
rw0=[rw12,rw13,rw14,rw23,rw24,rw34];
rpwall=[rpw12;rpw13;rpw14;rpw23;rpw24;rpw34];
pval=[pvalue12,pvalue13,pvalue14,pvalue23,pvalue24,pvalue34];
