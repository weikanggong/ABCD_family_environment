clear
close
clc

addpath(genpath('.../matlab_toolbox/CanlabCore-master'))
addpath(genpath('.../matlab_toolbox/MediationToolbox-master'))

load('family_data_for_mediation_19Mar2020.mat')
brain_id = id;

%% %%%%%%%%%%%%%% find common subjects with interested pheno
load freesurfer_data_for_region_analysis_3Aug2019

pheno_path = '.../';

instere_pheno = {'cognition'; 'depression'; 'psytotal'};
pheno_intere_id = 3;

intere_pheno = 'abcd_sscey01';
load([pheno_path, intere_pheno ,'.mat']);
eval(['pheno_8 = ' intere_pheno ';'])
pheno_8_bl = pheno_8(cellfun(@isempty,strfind(pheno_8.eventname,'baseline_year_1_arm_1'))==0,:);
pheno_8_id = pheno_8_bl.subjectkey;
pheno_8_fu = pheno_8(cellfun(@isempty,strfind(pheno_8.eventname,'baseline_year_1_arm_1'))==1,:);
pheno_8_fu_id = pheno_8_fu.subjectkey;

if pheno_intere_id==1
    %%%%%%%%% cognition
    intere_pheno = 'abcd_tbss01';
    load([pheno_path, intere_pheno ,'.mat']);
    eval(['pheno_9 = ' intere_pheno ';'])
    pheno_9_bl = pheno_9(cellfun(@isempty,strfind(pheno_9.eventname,'baseline_year_1_arm_1'))==0,:);
    pheno_9_id = pheno_9_bl.subjectkey;
elseif pheno_intere_id==2
    %%%%%%%%% depression
    intere_pheno = 'abcd_cbcls01';
    load([pheno_path, intere_pheno ,'.mat']);
    eval(['pheno_9 = ' intere_pheno ';'])
    pheno_9_bl = pheno_9(cellfun(@isempty,strfind(pheno_9.eventname,'baseline_year_1_arm_1'))==0,:);
    pheno_9_id = pheno_9_bl.subjectkey;
else
    %%%%%%%%% depression
    intere_pheno = 'abcd_cbcls01';
    load([pheno_path, intere_pheno ,'.mat']);
    eval(['pheno_9 = ' intere_pheno ';'])
    pheno_9_bl = pheno_9(cellfun(@isempty,strfind(pheno_9.eventname,'baseline_year_1_arm_1'))==0,:);
    pheno_9_id = pheno_9_bl.subjectkey;
end

final_id_new = intersect_multi({final_id; pheno_8_id; pheno_9_id; brain_id});

[~, index] = intersect(brain_id,final_id_new);
brain_data_final = area_para_monitor(index,:);

[~, index] = intersect(pheno_8_id,final_id_new);
pheno_8_final = pheno_8_bl(index,:);
%%%%%%%%%%
% envi_pheno_data = str2double(pheno_8_final.fes_y_ss_fc);
envi_pheno_data = str2double(pheno_8_final.pmq_y_ss_mean);
%%%%%%%%%%

[~, index] = intersect(pheno_9_id,final_id_new);
pheno_9_final = pheno_9_bl(index,:);

if pheno_intere_id==1
    intere_pheno_data = str2double(pheno_9_final.nihtbx_totalcomp_uncorrected);  %%%% cognition
elseif pheno_intere_id==2
    intere_pheno_data = str2double(pheno_9_final.cbcl_scr_dsm5_depress_r);   %%%% depression
else
    intere_pheno_data = str2double(pheno_9_final.cbcl_scr_syn_totprob_r);   %%%% total psychiatry
end

[~, index] = intersect(final_id, final_id_new);
cov_data = [age, sex, educ_data, income_data, bmi, race_1, race_2, race_3, puberty];
cov_data = cov_data(index,:);

site_data_1 = str2double(strrep(site_data_final(index),'site',''));
site_data_final = dummyvar(site_data_1);
site_data_final = site_data_final(:,2:end);

cov_data = [cov_data, site_data_final];

fam_stru = str2double(fam_stru(index));
[bb,~] = unique(fam_stru);

index_fam = nan(length(fam_stru),1);
for i=1:length(bb)
    index_fam_1 = find(fam_stru==bb(i));
    if length(index_fam_1)==1
        index_fam(index_fam_1,1) = 1;
    elseif length(index_fam_1)==2
        index_fam(index_fam_1(1),1) = 1;
        index_fam(index_fam_1(2),1) = 0;
    else
        index_fam(index_fam_1(1),1) = 1;
        index_fam(index_fam_1(2),1) = 0;
        index_fam(index_fam_1(3),1) = 0;
    end
end

%% %%%%%%%%%%%%%%%%
X = brain_data_final;
M = intere_pheno_data;
Y = envi_pheno_data;

index_nan = sum([X, Y, M, cov_data],2);
X_data = X(~isnan(index_nan));
M_data = M(~isnan(index_nan));
Y_data = Y(~isnan(index_nan));
cov_all_data = cov_data(~isnan(index_nan),:);

%%%%%%%%%%%%%%
b = regress(X_data,[cov_all_data ones(size(cov_all_data,1),1)]);
X_data_cov = X_data - (b(1:end-1)'*cov_all_data')';

b = regress(M_data,[cov_all_data ones(size(cov_all_data,1),1)]);
M_data_cov = M_data - (b(1:end-1)'*cov_all_data')';

b = regress(Y_data,[cov_all_data ones(size(cov_all_data,1),1)]);
Y_data_cov = Y_data - (b(1:end-1)'*cov_all_data')';

%%%%%%%%%%%%%%%

[paths_mean1, stats_mean1] = mediation(zscore(X_data_cov),zscore(Y_data_cov),zscore(M_data_cov),'boottop','bootsamples',10000,'doCIs');

[paths_mean2, stats_mean2] = mediation(zscore(X_data_cov),zscore(Y_data_cov),zscore(M_data_cov));












