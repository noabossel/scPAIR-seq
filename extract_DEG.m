
function [up_reg_DEGs,down_reg_DEGs]=extract_DEG(curr_mut)

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

if ~ismember(curr_mut,mut_list)
    error('Name of mutant is not valid; e.g., sifA, sifB, ssaV')
end

%generate list for each mutant (up/down)
sig_genes_per_mut_up=cell(n_mut,1);
sig_genes_per_mut_down=cell(n_mut,1);
all_genes_list=[];
for i=1:n_mut
    curr_mut_tf=ismember(mut_name,mut_list{i});
    curr_scores=maestro_scores(curr_mut_tf);
    curr_genes=gene_name_for_scores(curr_mut_tf);
    tf_in=curr_scores>=thresh_95;
    all_genes_list=[all_genes_list;curr_genes(tf_in)];
    curr_genes_mean_WT=mean_log2_wt(curr_mut_tf);
    curr_genes_mean_mut=mean_log2_mut(curr_mut_tf);
    curr_genes_fc=curr_genes_mean_mut-curr_genes_mean_WT;
    sig_genes_per_mut_up{i}=curr_genes(tf_in&curr_genes_fc>0);
    sig_genes_per_mut_down{i}=curr_genes(tf_in&curr_genes_fc<0);
end

%check for each gene how many times appear:
all_sig_genes_uniq=unique(all_genes_list);
n_genes=length(all_sig_genes_uniq);
tf_gene_up_in_mut=zeros(n_genes,n_mut);
tf_gene_down_in_mut=zeros(n_genes,n_mut);
for i=1:n_genes
    curr_gene=all_sig_genes_uniq{i};
    curr_gene_up_i=find(cellfun(@(x) sum(ismember(x,curr_gene)), sig_genes_per_mut_up));
    if ~isempty(curr_gene_up_i)
        tf_gene_up_in_mut(i,curr_gene_up_i)=1;
    end
    curr_gene_down_i=find(cellfun(@(x) sum(ismember(x,curr_gene)), sig_genes_per_mut_down));
    if ~isempty(curr_gene_down_i)
        tf_gene_down_in_mut(i,curr_gene_down_i)=1;
    end
end
n_time_each_gene_appear_up_down=[sum(tf_gene_up_in_mut,2),sum(tf_gene_down_in_mut,2)];
%genes is 10 mutants or more:
n_shared=10;
genes_up_shared=all_sig_genes_uniq(n_time_each_gene_appear_up_down(:,1)>=n_shared);
genes_down_shared=all_sig_genes_uniq(n_time_each_gene_appear_up_down(:,2)>=n_shared);

%specifc list for mutant:
i=strmatch(curr_mut,mut_list);
genes_up_all=sig_genes_per_mut_up{i};
up_reg_DEGs=setdiff(genes_up_all,genes_up_shared);
genes_down_all=sig_genes_per_mut_down{i};
down_reg_DEGs=setdiff(genes_down_all,genes_down_shared);

writecell(up_reg_DEGs,['output_files/',curr_mut,'_up_DEG.txt'],'Delimiter','tab')
writecell(down_reg_DEGs,['output_files/',curr_mut,'_down_DEG.txt'],'Delimiter','tab')

