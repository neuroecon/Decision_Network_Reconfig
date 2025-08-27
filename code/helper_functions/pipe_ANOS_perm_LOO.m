% pipeline of anosperm
function [rw0,rpwall,pval]=pipe_ANOS_perm_LOO(datadir,perm)

% rw0: actual stat
% rpall: all the stats after permutation
% pval: p-value

[rw12,rpw12,pvalue12]=ANOS_perm(datadir,1,2,perm);
[rw13,rpw13,pvalue13]=ANOS_perm(datadir,1,3,perm);
[rw14,rpw14,pvalue14]=ANOS_perm(datadir,1,4,perm);
[rw23,rpw23,pvalue23]=ANOS_perm(datadir,2,3,perm);
[rw24,rpw24,pvalue24]=ANOS_perm(datadir,2,4,perm);
[rw34,rpw34,pvalue34]=ANOS_perm(datadir,3,4,perm);

rw0=[rw12,rw13,rw14,rw23,rw24,rw34];
rpwall=[rpw12;rpw13;rpw14;rpw23;rpw24;rpw34];
pval=[pvalue12,pvalue13,pvalue14,pvalue23,pvalue24,pvalue34];
