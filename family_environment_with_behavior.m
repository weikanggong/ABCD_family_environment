clear
close
clc

%%{
%% %%%%%%%%%%%%%% load all data
load freesurfer_data_for_region_analysis_3Aug2019

intere_pheno = 'abcd_sscey01';
load([pheno_path, intere_pheno ,'.mat']);
eval(['pheno_9 = ' intere_pheno ';'])
pheno_9_bl = pheno_9(cellfun(@isempty,strfind(pheno_9.eventname,'baseline_year_1_arm_1'))==0,:);
pheno_9_id = pheno_9_bl.subjectkey;
pheno_9_fu = pheno_9(cellfun(@isempty,strfind(pheno_9.eventname,'baseline_year_1_arm_1'))==1,:);
pheno_9_fu_id = pheno_9_fu.subjectkey;

final_id_new = intersect_multi({final_id; pheno_9_id});

[~, index] = intersect(pheno_9_id,final_id_new);
pheno_9_final = pheno_9_bl(index,:);
%%%%%%%%%%
% select_pheno = str2double(pheno_9_final.fes_y_ss_fc);
select_pheno = str2double(pheno_9_final.pmq_y_ss_mean);
%%%%%%%%%%

[~, index] = intersect(final_id, final_id_new);
cov_data = [age, sex, educ_data, income_data, bmi, race_1, race_2, race_3, puberty];
cov_data = cov_data(index,:);
brain_data_final = brain_data_final(index,:);

site_id = str2double(strrep(site_data_final(index),'site',''));
fam_stru = str2double(fam_stru(index));

%% %%%%%%%%%%%%%
delete(gcp('nocreate'))
parpool(20)

dir_path = dir([pheno_path, '*.mat']);

count = 1;
for i=1:length(dir_path)
    i
    file_name = dir_path(i).name;
    if isempty(dir(['..../pmq_y_ss_mean_pheno_result/' file_name(1:end-4) '_result.mat']))
        
        load([pheno_path, file_name])
        eval(['other_pheno = ' file_name(1:end-4) ';'])
        
        phen_leng(i,1) = size(other_pheno,2);
        pheno_name_index{i,1} = dir_path(i).name(1:end-4);
        
        if sum(strcmpi(other_pheno.Properties.VariableNames,'eventname'))~=0
            other_pheno = other_pheno(cellfun(@isempty,strfind(other_pheno.eventname,'baseline_year_1_arm_1'))==0,:);
        end
        
        if size(other_pheno,1)>1000 && size(other_pheno,2)<1000 ...
                other_pheno_id = other_pheno.subjectkey;
            other_pheno_name = other_pheno.Properties.VariableNames';
            
            final_id_1 = intersect(other_pheno_id,final_id_new);
            [~, index_pheno] = intersect(final_id_new,final_id_1);
            select_pheno_final = select_pheno(index_pheno,:);
            cov_data_final = cov_data(index_pheno,:);
            site_id_final = site_id(index_pheno,:);
            fam_stru_final = fam_stru(index_pheno,:);
            
            [~, index_pheno] = intersect(other_pheno_id,final_id_1);
            other_pheno_final = other_pheno(index_pheno,:);
            
            tval = nan(size(other_pheno_final,2),1);
            rval = nan(size(other_pheno_final,2),1);
            pval = nan(size(other_pheno_final,2),1);
            parfor j=9:size(other_pheno_final,2)
                if size(other_pheno_final(:,j),1)>2000 && sum(isnan(str2double(table2array(other_pheno_final(:,j)))))<9000 ...
                        && sum(unique(str2double(table2array(other_pheno_final(:,j))))>=0)>1
                    try
                        index = sum([select_pheno_final,cov_data_final, str2double(table2array(other_pheno_final(:,j)))],2);
                        data_lme = array2table([select_pheno_final,cov_data_final]);
                        data_lme.Properties.VariableNames = {'sleep','age', 'sex', 'educ', 'income', 'bmi', 'race_1', 'race_2', 'race_3','pertu'};
                        data_lme.site = site_id_final;
                        data_lme.fam_stru = fam_stru_final;
                        data_lme.smri_data1 = str2double(table2array(other_pheno_final(:,j)));
                        data_lme(isnan(index),:) = [];
                        
                        lme = fitlme(data_lme,'smri_data1~sleep+age+sex+educ+income+bmi+race_1+race_2+race_3+pertu+(1|site:fam_stru)');
                        
                        num_sub = sum(~isnan(data_lme.smri_data1));
                        
                        tval(j,1) = table2array(lme.Coefficients(2,4));
                        rval(j,1) = tval(j,1)/sqrt(num_sub-2+tval(j,1)*tval(j,1));
                        pval(j,1) = table2array(lme.Coefficients(2,6));
                    catch
                    end
                end
            end
            
            eval(['save ..../pmq_y_ss_mean_pheno_result/' file_name(1:end-4) '_result.mat tval rval pval other_pheno_name'])
            
%             cval_all{i,1} = cval;
%             pval_all{i,1} = pval;
%             pheno_name_all_full{i,1} = other_pheno_name;
%             pheno_name_all{i,1} = dir_path(i).name(1:end-4);
        end
        
        eval(['clear ' file_name(1:end-4)])
    end
end

% save all_pheno_vs_sleep_duration_2Aug2019.mat cval_all pval_all pheno_name_all pheno_name_all_full





