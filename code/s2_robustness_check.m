% robustness check: changing density and FC construction methods
wd = pwd;
addpath([wd,'/helper_functions']);
color_map = {'#808080','#f7c334','#ed98b3','#ba93db'}; % for 4 conditions
color_map1 = {[128, 128, 128]/256, [247, 195, 52]/256,[237, 152, 179]/256,[187, 147, 219]/256};
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

%% Partition similarity between different param choice and default option
resultDir = '../results/community_detection/';

% calculate similarity within an individual iteration
subList=[702 705 708 711 718 719 720 722 725 728 729 730 733 735 ...
736 737 740 744 745 748 750 754 755 758 759 760 761];
MIn_sub = struct();
for sub=1:length(subList)
    foldername = ['Sub',num2str(subList(sub))];
    Sall_005 = load([resultDir,'/0.05,1,0.1_concat/',foldername,'/S100.mat']);
    Sall_010 = load([resultDir,'/0.1,1,0.1_concat/',foldername,'/S100.mat']);
    Sall_020 = load([resultDir,'/0.2,1,0.1_concat/',foldername,'/S100.mat']);

    Sall_005_avg = load([resultDir,'/0.05,1,0.1_avg/',foldername,'/S100.mat']);
    Sall_010_avg = load([resultDir,'/0.1,1,0.1_avg/',foldername,'/S100.mat']);
    Sall_020_avg = load([resultDir,'/0.2,1,0.1_avg/',foldername,'/S100.mat']);

    for rep=1:100
        for j=1:4
            [~, MIn_sub.den005(sub,rep,j)] = partition_distance(Sall_005.Sallps(:,rep,j),Sall_010.Sallps(:,rep,j));
            [~, MIn_sub.den020(sub,rep,j)] = partition_distance(Sall_020.Sallps(:,rep,j),Sall_010.Sallps(:,rep,j));
            [~, MIn_sub.den005_avg(sub,rep,j)] = partition_distance(Sall_005_avg.Sallps(:,rep,j),Sall_010.Sallps(:,rep,j));
            [~, MIn_sub.den020_avg(sub,rep,j)] = partition_distance(Sall_020_avg.Sallps(:,rep,j),Sall_010.Sallps(:,rep,j));
            [~, MIn_sub.den010_avg(sub,rep,j)] = partition_distance(Sall_010_avg.Sallps(:,rep,j),Sall_010.Sallps(:,rep,j));
            [~, MIn_sub.den005_ac(sub,rep,j)] = partition_distance(Sall_005_avg.Sallps(:,rep,j),Sall_005.Sallps(:,rep,j));
            [~, MIn_sub.den020_ac(sub,rep,j)] = partition_distance(Sall_020_avg.Sallps(:,rep,j),Sall_020.Sallps(:,rep,j));

        end
    end
end
mean_sub_MIn = struct();
mean_sub_MIn.den005 = reshape(mean(MIn_sub.den005,2),27,4);
mean_sub_MIn.den020 = reshape(mean(MIn_sub.den020,2),27,4);
mean_sub_MIn.den005_avg = reshape(mean(MIn_sub.den005_avg,2),27,4);
mean_sub_MIn.den020_avg = reshape(mean(MIn_sub.den020_avg,2),27,4);
mean_sub_MIn.den010_avg = reshape(mean(MIn_sub.den010_avg,2),27,4);
mean_sub_MIn.den005_ac = reshape(mean(MIn_sub.den005_ac,2),27,4);
mean_sub_MIn.den020_ac = reshape(mean(MIn_sub.den020_ac,2),27,4);

mean_sub_MIn_all = [mean(mean_sub_MIn.den005,1);mean(mean_sub_MIn.den020,1);...
    mean(mean_sub_MIn.den005_avg,1);mean(mean_sub_MIn.den020_avg,1);...
    mean(mean_sub_MIn.den010_avg,1);mean(mean_sub_MIn.den005_ac,1);...
    mean(mean_sub_MIn.den020_ac,1)];

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


%% Comparison between concatenation and average
% compare the FC matrix
subList=[702 705 708 711 718 719 720 722 725 728 729 730 733 735 ...
736 737 740 744 745 748 750 754 755 758 759 760 761];

dataDir='../data/';
condlist = {'SF','EMC','IMC'};

FC_concat_all = {};
FC_avg_all = {};
E_concat = [];
E_avg = [];
lambda_concat = [];
lambda_avg = [];
Q_concat = [];
Q_avg = [];

for i=1:3
    cond = condlist{i};
    for sub = 1:length(subList)
        FC_mat_concat = load([dataDir,cond,'/FCmat_concat/zSub',num2str(subList(sub)),'.txt']);
        load([dataDir,cond,'/FCmat_avg/zSub',num2str(subList(sub)),'.mat']);
        % save all FC
        FC_concat_all{i}(:,:,sub) = FC_mat_concat;
        FC_avg_all{i}(:,:,sub) = FC_mat_avg;     
        % 10% network
        W_concat = threshold_proportional(FC_mat_concat, 0.1);
        W_avg = threshold_proportional(FC_mat_avg, 0.1);
        % calculate network efficiency
        E_concat(i,sub) = efficiency_wei(W_concat, 0);
        E_avg(i,sub) = efficiency_wei(W_avg,0);
        % calculate clustering coefficient
        lambda_concat(i,sub)= charpath(W_concat);
        lambda_avg(i,sub) = charpath(W_avg);
        % modularity
        [~,Q_concat(i,sub)]=modularity_und(W_concat,1);
        [~,Q_avg(i,sub)] = modularity_und(W_avg,1);
                
    end
end

FC_concat_SF = mean(FC_concat_all{1},3);
FC_avg_SF = mean(FC_avg_all{1},3);
FC_diff_SF = FC_concat_SF - FC_avg_SF;

% global efficiency
figure;
hold on
for i= [2,1,3]
    scatter(E_concat(i,:),E_avg(i,:),40,color_map1{i+1},'filled','MarkerFaceAlpha',0.7);
end
rl = refline([1,0]);
rl.Color = 'k';
title('Network Efficiency','Fontsize',15);
xlabel('Concatenation','Fontsize',13);
ylabel('Average','Fontsize',13);
leg = legend('EMC','SF','IMC');
set(leg,'box' , 'off')
set(gca,'fontname','Avenir') 
set(gcf,'position',[10,10,350,300])  
xlim([0.1,0.22]);
ylim([0.1,0.22]);
% print correlation
[r,p] = corr(E_concat(:),E_avg(:));
disp(['R = ',num2str(round(r,3)),', p = ',num2str(p)]);


% characteristic path length
figure;
hold on
for i=[2,1,3]
    scatter(lambda_concat(i,:),lambda_avg(i,:),40,color_map1{i+1},'filled','MarkerFaceAlpha',0.7);
end
rl = refline([1,0]);
rl.Color = 'k';
title('Characteristic Path Length','Fontsize',15);
xlabel('Concatenation','Fontsize',13);
ylabel('Average','Fontsize',13);
leg = legend('EMC','SF','IMC');
set(leg,'box' , 'off')
set(gca,'fontname','Avenir') 
set(gcf,'position',[10,10,350,300])  
% print correlation
[r,p] = corr(lambda_concat(:),lambda_avg(:));
disp(['R = ',num2str(round(r,3)),', p = ',num2str(p)]);


% modularity
figure;
hold on
for i=[2,1,3]
    scatter(Q_concat(i,:),Q_avg(i,:),40,color_map1{i+1},'filled','MarkerFaceAlpha',0.7);
end
rl = refline([1,0]);
rl.Color = 'k';
title('Modularity','Fontsize',15);
xlabel('Concatenation','Fontsize',13);
ylabel('Average','Fontsize',13);
leg = legend('EMC','SF','IMC');
set(leg,'box' , 'off')
set(gca,'fontname','Avenir') 
set(gcf,'position',[10,10,350,300])  
% print correlation
[r,p] = corr(Q_concat(:),Q_avg(:));
disp(['R = ',num2str(round(r,3)),', p = ',num2str(p)]);



% nodal degree
figure;
hold on
for i = [2,1,3]
    deg_concat = degrees_und(threshold_proportional(mean(FC_concat_all{i},3), 0.1));
    deg_avg = degrees_und(threshold_proportional( mean(FC_avg_all{i},3), 0.1));
    scatter(deg_concat,deg_avg,30,color_map1{i+1},'filled','MarkerFaceAlpha',0.5);
end
rl = refline([1,0]);
rl.Color = 'k';
title('Nodal Degree','Fontsize',15);
xlabel('Concatenation','Fontsize',13);
ylabel('Average','Fontsize',13);
leg = legend('EMC','SF','IMC');
set(leg,'box' , 'off')
set(gca,'fontname','Avenir') 
set(gcf,'position',[10,10,350,300])  
% print correlation
[r,p] = corr(deg_concat(:),deg_avg(:));
disp(['R = ',num2str(round(r,3)),', p = ',num2str(p)]);


% nodal clustering coefficient
figure;
hold on
for i = [2,1,3]
    cc_concat = clustering_coef_wu(threshold_proportional(mean(FC_concat_all{i},3), 0.1));
    cc_avg = clustering_coef_wu(threshold_proportional( mean(FC_avg_all{i},3), 0.1));
    scatter(cc_concat,cc_avg,30,color_map1{i+1},'filled','MarkerFaceAlpha',0.5);
end
rl = refline([1,0]);
rl.Color = 'k';
title('Nodal Clustering Coefficient','Fontsize',15);
xlabel('Concatenation','Fontsize',13);
ylabel('Average','Fontsize',13);
leg = legend('EMC','SF','IMC');
set(leg,'box' , 'off')
set(gca,'fontname','Avenir') 
set(gcf,'position',[10,10,350,300])  
xlim([0,0.6])
% print correlation
[r,p] = corr(cc_concat(:),cc_avg(:));
disp(['R = ',num2str(round(r,3)),', p = ',num2str(p)]);



%% LOO global level between condition partition smilarity test
% 1. community detection
subList=[702 705 708 711 718 719 720 722 725 728 729 730 733 735 ...
736 737 740 744 745 748 750 754 755 758 759 760 761];

dataDir='../data/';
resultDir = '../results/';

parfor (i=[10:21],8)
    subList_new = subList;
    subList_new(i) = [];
    sub_exc = subList(i);
    community_detection_global_loo(dataDir,resultDir,'concat',0.1,1,0.1,subList_new,num2str(sub_exc));
end

% 2. calculate Variation of Information (VIn) as a quantitative measure of
% network similairty
subList=[702 705 708 711 718 719 720 722 725 728 729 730 733 735 ...
736 737 740 744 745 748 750 754 755 758 759 760 761];
datadir = '../results/community_detection/LOO/0.1,1,0.1_concatno_';
VIn_sub = {};
VIn = [];
for sub=1:length(subList)
    load([datadir,num2str(subList(sub)),'/S27.mat']);
    for s=1:size(Sall,2)
        for i=1:3
            for j=i+1:4
                [VIn_sub{sub}(i,j,s),~] = partition_distance(Sall(:,s,i), Sall(:,s,j));
            end
        end
    end
    VIn(:,:,sub) = mean(VIn_sub{sub},3);
end
VIn_mean = mean(VIn, 3); % mean
VIn_std = std(VIn,[],3); % STD

% 3. perform Multi-Response Permutation Procedure to test the statistical
% significance of the differences
perm = 5000; % number of permutations
pval = [];
parfor (sub = 1:length(subList),12)
    subname = num2str(subList(sub));
    [~,~,pval(sub,:)]=pipe_ANOS_perm_LOO([datadir,subname],perm);
end
    

%% LOO local level: between-condition community similarity analysis
% 1. community detection 
subList=[702 705 708 711 718 719 720 722 725 728 729 730 733 735 ...
736 737 740 744 745 748 750 754 755 758 759 760 761];

dataDir='../data/';
resultDir = '../results/';
load([resultDir,'community_detection/0.1,1,0.1_concat/Subnetwork/roilist.mat']);

parfor (i=1:27,8)
    subList_new = subList;
    subList_new(i) = [];
    sub_exc = subList(i);
    community_detection_local_loo(dataDir,resultDir,roilist,'concat',0.1,1,0.1,subList_new,num2str(sub_exc));
end

% 2. calculate Variation of Information (VIn) as a quantitative measure of
% network similairty
subList=[702 705 708 711 718 719 720 722 725 728 729 730 733 735 ...
736 737 740 744 745 748 750 754 755 758 759 760 761];
datadir = '../results/community_detection/LOO/0.1,1,0.1_concatno_';
VIn_sub = {};
VIn = [];
for sub=1:length(subList)
    load([datadir,num2str(subList(sub)),'/Subnetwork/S27.mat']);
    for s=1:size(Sall,2)
        for i=1:3
            for j=i+1:4
                [VIn_sub{sub}(i,j,s),~] = partition_distance(Sall(:,s,i), Sall(:,s,j));
            end
        end
    end
    VIn(:,:,sub) = mean(VIn_sub{sub},3);
end
VIn_mean = mean(VIn, 3); % mean
VIn_std = std(VIn,[],3); % STD


% 3. perform Multi-Response Permutation Procedure to test the statistical
% significance of the differences
perm = 5000; % number of permutations
pval = [];
parfor (sub = 1:length(subList),8)
    subname = num2str(subList(sub));
    [~,~,pval(sub,:)]=pipe_ANOS_perm_LOO([datadir,subname,'/Subnetwork'],perm);
end

%% change gamma
dataDir='../data/';
resultDir = '../results/';
density = 0.1;
omega = 0.1;
gamma_list = [0.25,0.5,0.75,1.25,1.5,1.75]; 
subList=[702 705 708 711 718 719 720 722 725 728 729 730 733 735 ...
736 737 740 744 745 748 750 754 755 758 759 760 761];

% global level community detection
parfor (g = 1:6,6)
    gamma = gamma_list(g);
    community_detection_global(dataDir,resultDir,'concat',density,gamma,omega);
end

% calculate similarity within an individual iteration
MIn_sub = {};
for sub=1:length(subList)
    for g = 1:6
        gamma = gamma_list(g);
        params = ['0.1,',num2str(gamma),',0.1_concat/'];
        foldername = ['Sub',num2str(subList(sub))];
        Sall_gamma = load([resultDir,'community_detection/',params,foldername,'/S100.mat']);
        Sall_template = load([resultDir,'community_detection/0.1,1,0.1_concat/',foldername,'/S100.mat']);
        for rep=1:100
            for j=1:4
                [~, MIn_sub{g}(sub,rep,j)] = partition_distance(Sall_gamma.Sallps(:,rep,j),Sall_template.Sallps(:,rep,j));
            end
        end
    end
end
mean_sub_MIn = {};
mean_sub_MIn_all = [];
E = [];
for g = 1:6
    mean_sub_MIn{g}= reshape(mean(MIn_sub{g},2),27,4);
    mean_sub_MIn_all(g,:) = mean(mean_sub_MIn{g},1);
    E(g,:) = std(mean_sub_MIn{g},[],1)./sqrt(27);
end

% plot the partition similarity results
figure; hold on;
for k = [1,3,2,4]
    errorbar(gamma_list, mean_sub_MIn_all(:,k), E(:,k), '-', ...
        'LineWidth', 1.2,'Color',color_map{k},'CapSize',0);  % one line with error bars
end
hold off;
leg = legend('Task condition');
set(leg,'box' , 'off')
set(gca,'fontname','Avenir') 
set(gcf,'position',[10,10,500,350])  
legend('Rest','EMC','SF','IMC');
xlabel('Structural Resolution \gamma','FontSize',13);
ylabel('MIn with Default Partition','FontSize',13);
ylim([0.2,1]);

% local level community detection
load([resultDir,'community_detection/0.1,1,0.1_concat/Subnetwork/roilist.mat']);
parfor (g=1:2,2)
    gamma = gamma_list(g);
    community_detection_local(dataDir,resultDir,roilist, 'concat',density,gamma,omega);
end

% local level similarity
MIn_sub = {};
for sub=1:length(subList)
    for g = 1:6
        gamma = gamma_list(g);
        params = ['0.1,',num2str(gamma),',0.1_concat/'];
        foldername = ['Sub',num2str(subList(sub))];
        Sall_gamma = load([resultDir,'community_detection/',params,'Subnetwork/',foldername,'/localS100.mat']);
        Sall_template = load([resultDir,'community_detection/0.1,1,0.1_concat/Subnetwork/',foldername,'/118S100.mat']);
        for rep=1:100
            for j=1:4
                [~, MIn_sub{g}(sub,rep,j)] = partition_distance(Sall_gamma.Sallps(:,rep,j),Sall_template.Sallps(:,rep,j));
            end
        end
    end
end
mean_sub_MIn = {};
mean_sub_MIn_all = [];
E = [];
for g = 1:6
    mean_sub_MIn{g}= reshape(mean(MIn_sub{g},2),27,4);
    mean_sub_MIn_all(g,:) = mean(mean_sub_MIn{g},1);
    E(g,:) = std(mean_sub_MIn{g},[],1)./sqrt(27);
end

% plot the partition similarity results
figure; hold on;
for k = [1,3,2,4]
    errorbar(gamma_list, mean_sub_MIn_all(:,k), E(:,k), '-', ...
        'LineWidth', 1.2,'Color',color_map{k},'CapSize',0);  % one line with error bars
end
hold off;
leg = legend('Task condition');
set(leg,'box' , 'off')
set(gca,'fontname','Avenir') 
set(gcf,'position',[10,10,500,350])  
legend('Rest','EMC','SF','IMC');
xlabel('Structural Resolution \gamma','FontSize',13);
ylabel('MIn with Default Partition','FontSize',13);
ylim([0.2,1]);

%% change omega
color_map = {'#808080','#f7c334','#ed98b3','#ba93db'};
dataDir='../data/';
resultDir = '../results/';
density = 0.1;
gamma = 1;
omega_list = [0.05,0.15, 0.2, 0.3,0.5, 0.6,0.9,1.2,1.5];
subList=[702 705 708 711 718 719 720 722 725 728 729 730 733 735 ...
736 737 740 744 745 748 750 754 755 758 759 760 761];

% global level community detection
parfor (o = 1:5,5)
    omega = omega_list(o);
    community_detection_global(dataDir,resultDir,'concat',density,gamma,omega);
end

% calculate similarity within an individual iteration
MIn_sub = {};
for sub=1:length(subList)
    for o = 1:length(omega_list)
        omega = omega_list(o);
        params = ['0.1,1,',num2str(omega),'_concat/'];
        foldername = ['Sub',num2str(subList(sub))];
        Sall_gamma = load([resultDir,'community_detection/',params,foldername,'/S100.mat']);
        Sall_template = load([resultDir,'community_detection/0.1,1,0.1_concat/',foldername,'/S100.mat']);
        for rep=1:100
            for j=1:4
                [~, MIn_sub{o}(sub,rep,j)] = partition_distance(Sall_gamma.Sallps(:,rep,j),Sall_template.Sallps(:,rep,j));
            end
        end
    end
end
mean_sub_MIn = {};
mean_sub_MIn_all = [];
E = [];
for o = 1:length(omega_list)
    mean_sub_MIn{o}= reshape(mean(MIn_sub{o},2),27,4);
    mean_sub_MIn_all(o,:) = mean(mean_sub_MIn{o},1);
    E(o,:) = std(mean_sub_MIn{o},[],1)./sqrt(27);
end

% plot the partition similarity results
figure; hold on;
for k = [1,3,2,4]
    errorbar(omega_list, mean_sub_MIn_all(:,k), E(:,k), '-', ...
        'LineWidth', 1.2,'Color',color_map{k},'CapSize',0);  % one line with error bars
end
hold off;
leg = legend('Task condition');
set(leg,'box' , 'off')
set(gca,'fontname','Avenir') 
set(gcf,'position',[10,10,500,350])  
legend('Rest','EMC','SF','IMC');
xlabel('Temporal Resolution \omega','FontSize',13);
ylabel('MIn with Default Partition','FontSize',13);
ylim([0.2,1]);
xlim([0,1.6]);


% local level community detection
load([resultDir,'community_detection/0.1,1,0.1_concat/Subnetwork/roilist.mat']);
parfor (o = 1:5,5)
    omega = omega_list(o);
    community_detection_local(dataDir,resultDir,roilist, 'concat',density,gamma,omega);
end

% calculate similarity within an individual iteration
MIn_sub = {};
for sub=1:length(subList)
    for o = 1:length(omega_list)
        omega = omega_list(o);
        params = ['0.1,1,',num2str(omega),'_concat/'];
        foldername = ['Sub',num2str(subList(sub))];
        Sall_omega = load([resultDir,'community_detection/',params,'Subnetwork/',foldername,'/localS100.mat']);
        Sall_template = load([resultDir,'community_detection/0.1,1,0.1_concat/Subnetwork/',foldername,'/118S100.mat']);
        for rep=1:100
            for j=1:4
                [~, MIn_sub{o}(sub,rep,j)] = partition_distance(Sall_omega.Sallps(:,rep,j),Sall_template.Sallps(:,rep,j));
            end
        end
    end
end
mean_sub_MIn = {};
mean_sub_MIn_all = [];
E = [];
for o = 1:length(omega_list)
    mean_sub_MIn{o}= reshape(mean(MIn_sub{o},2),27,4);
    mean_sub_MIn_all(o,:) = mean(mean_sub_MIn{o},1);
    E(o,:) = std(mean_sub_MIn{o},[],1)./sqrt(27);
end

% plot the partition similarity results
figure; hold on;
for k = [1,3,2,4]
    errorbar(omega_list, mean_sub_MIn_all(:,k), E(:,k), '-', ...
        'LineWidth', 1.2,'Color',color_map{k},'CapSize',0);  % one line with error bars
end
hold off;
leg = legend('Task condition');
set(leg,'box' , 'off')
set(gca,'fontname','Avenir') 
set(gcf,'position',[10,10,500,350])  
legend('Rest','EMC','SF','IMC');
xlabel('Temporal Resolution \omega','FontSize',13);
ylabel('MIn with Default Partition','FontSize',13);
ylim([0.2,1]);
xlim([0,1.6]);
