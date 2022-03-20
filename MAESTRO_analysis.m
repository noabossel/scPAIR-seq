
clear all


close all
clc

%load real MAESTRO scores:
data=importdata('maestro_scores_real_data.txt');
mut_name=data.textdata(2:end,3);
gene_name_for_scores=data.textdata(2:end,4);
maestro_scores=data.data(:,1);
mean_log2_wt=data.data(:,7);
mean_log2_mut=data.data(:,2);
diff_wt_mut=mean_log2_mut-mean_log2_wt;

thresh_95=prctile(maestro_scores,95);

mut_list=unique(mut_name);
n_mut=length(mut_list);

% %upload random data 1000 permutations:
data=importdata('maestro_scores_rand.txt');
mut_name_rand=data.textdata(2:end,2);
gene_name_for_maestro_rand=data.textdata(2:end,3);
maestro_scores_rand=data.data;

thresh_95_rand=prctile(maestro_scores_rand,95);

% find index of the start of each permutation
tmp_i=strmatch(mut_name_rand{1},mut_name_rand); %The first mutant appears in the begining of each permutation
diff_stee_i=diff(tmp_i);
i_new_perm=find(diff_stee_i>10)+1;
perm_i_indexes=[[1;tmp_i(i_new_perm)],[tmp_i(i_new_perm)-1;length(mut_name_rand)]];


n_perm=1000;
rand_maestro_scores_by_perm_by_mut=cell(n_perm,n_mut);
for i=1:n_perm
    curr_mut_name_rand=mut_name_rand(perm_i_indexes(i,1):perm_i_indexes(i,2));
    curr_maestro_data_rand=maestro_scores_rand(perm_i_indexes(i,1):perm_i_indexes(i,2));
    for j=1:n_mut
        curr_rand_mut_tf=ismember(curr_mut_name_rand,mut_list{j});
        curr_rand_maestro=curr_maestro_data_rand(curr_rand_mut_tf);
        rand_maestro_scores_by_perm_by_mut{i,j}=curr_rand_maestro;
    end
end

%calc p-values real vs. rand null dist. 
p_val_all=ones(n_mut,1);
diff_maestro_scores=zeros(n_mut,1);
for i=1:n_mut
    %real:
    curr_mut_tf=ismember(mut_name,mut_list{i});
    curr_maestro_mut=maestro_scores(curr_mut_tf);
    mean_maestro_real=mean(curr_maestro_mut);
    %rand all
    curr_all_rand_scores=[];
    for pi=1:n_perm
        curr_all_rand_scores=[curr_all_rand_scores;rand_maestro_scores_by_perm_by_mut{pi,i}];
        mean_maestro_rand=mean(curr_all_rand_scores);
    end
    [~,p_val_all(i)]=ttest2(curr_maestro_mut,curr_all_rand_scores);
    diff_maestro_scores(i)=mean_maestro_real-mean_maestro_rand;
end

q_val_all=mafdr(p_val_all);
sig_mut_tf=q_val_all<=0.01&diff_maestro_scores>0;

%summarize significant mutants:
sig_mutants=mut_list(sig_mut_tf);
p_val=p_val_all(sig_mut_tf);
q_val=q_val_all(sig_mut_tf);
out=table(sig_mutants,p_val,q_val);
writetable(out,'output_files/significant_mutants.txt');

