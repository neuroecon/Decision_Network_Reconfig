function [BetweenFCmean,BetweenFCsum,FCall]=BMstrength(FCmat,partition, mod1, mod2)
% between network connectivity
% FCmat: n*n FC matrix
% partition: n*1 vector
% mod1: visual vs task-evoked (can be either number or vector)
% mod2: task-evoked

FCall=[];

for i=1:size(FCmat,1)
    for j=1:size(FCmat,1)
        if ismember(partition(i), mod1) && ismember(partition(j),mod2)
            if FCmat(i,j)~=0
                FCall(end+1)=FCmat(i,j);  % remove 0 connections
            end
        end
    end
end

BetweenFCmean=mean(FCall)';
BetweenFCsum=sum(FCall)';

FCall=FCall';

