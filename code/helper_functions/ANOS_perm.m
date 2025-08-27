function  [rw,rpw,pvalue]=ANOS_perm(datadir,cond1,cond2,perm)
% permutation by anoism using a delta method
% condition: 1:rest 2:SF 3:EMC 4:IMC

load([datadir,'/S27.mat']);

n=size(Sall,2);
rwv1=[];
rwv2=[];
for j=1:n
    for k=j+1:n
        [rwv1(end+1),MIn]=partition_distance(Sall(:,j,cond1),Sall(:,k,cond1));
        [rwv2(end+1),MIn]=partition_distance(Sall(:,j,cond2),Sall(:,k,cond2));
    end
end
rw=mean([rwv1 rwv2]);

% permutation
subList = 1:n;
for q=1:perm
    rpv=[];
    % choose half of the subjects
    p = randperm(n,round(n/2));
    permlist=subList(p);
    temp1=Sall(:,:,cond1);
    temp2=Sall(:,:,cond2);
    for k=1:length(subList)
        if sum(ismember(p,k))
            temp1(:,k)=Sall(:,k,cond2);
            temp2(:,k)=Sall(:,k,cond1);
        end
    end
    
rpwv1=[];
rpwv2=[];
for j=1:n
    for k=j+1:n
        [rpwv1(end+1),MIn]=partition_distance(temp1(:,j),temp1(:,k));
        [rpwv2(end+1),MIn]=partition_distance(temp2(:,j),temp2(:,k));
    end
end
rpw(q)=mean([rpwv1 rpwv2]);
pvalue=1-findp(rpw,rw);
end

