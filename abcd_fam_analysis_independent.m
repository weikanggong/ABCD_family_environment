clear all

%% freesurfer
area_data=readNPY('.../abcd_area.npy');
fs_ids=readtable('.../abcd_fs_id.txt','ReadVariableNames',false,'Delimiter','\t');
fs_ids=fs_ids.Var1;

fs_ids1=fs_ids;
for i=1:10888
    f=fs_ids{i};
    fs_ids1{i}=strrep(f,'sub-NDAR','NDAR_');
end

fs_ids2=fs_ids1(area_data(:,1)~=0);
area_data=area_data(area_data(:,1)~=0,:);

%% select 1 child in each family
aa1=readtable('.../acspsw03.txt','delimiter','tab');
ids3=aa1.subjectkey(2:end);
fam_stru=str2double(aa1.rel_family_id(2:end));
fam_id=unique(fam_stru);

subj_id={};
for i=1:length(fam_id)
    
    indx=find(fam_stru==fam_id(i));
    subj_id{i,1}=ids3{indx(1)};    
end

[inter11,ia11,ib11]=intersect(fs_ids2,subj_id);

area_data1=area_data(ia11,:);

%% sepatate glm
cc=load('.../freesurfer_data_for_region_analysis.mat');
%cov_all = [age, sex, educ_data, income_data, bmi, race_1, race_2, race_3, site_data_final, birth_weight, early_birth, puberty];

aa11=readtable('..../fhxp201.txt','delimiter','tab');
ids2=aa11.subjectkey(2:end);
score01=aa11.fam_history_q11a_professional(2:end);
score02=aa11.fam_history_q11d_professional(2:end);
fam_hist=[score01,score02];

aa=readtable('.../abcd_sscey01.txt','delimiter','tab');
ids1=aa.subjectkey(2:end);
score1=(aa.fes_y_ss_fc(2:end));
score2=(aa.pmq_y_ss_mean(2:end));

inter=intersect_multi({ids1,inter11,cc.final_id,ids2});
[inter1,ia1]=intersect(ids1,inter);
[inter2,ia2]=intersect(inter11,inter);
[inter3,ia3]=intersect(cc.final_id,inter);
[inter4,ia4]=intersect(ids2,inter);

score11=str2double(score1(ia1,:));
score22=str2double(score2(ia1,:));

indx1=~isnan(score11) & ~isnan(score11);
indx2=~isnan(score11) & ~isnan(score11);

area_data1=area_data1(ia2,:);

cov_data=[g_impute_withmean(cc.cov_all(ia3,[1:28,31])),fam_hist1];

fam_hist1=double(str2double(fam_hist(ia4,:))==1);

Ts1=BWAS_Tregression([score11(indx1,:),cov_data(indx1,:)],area_data1(indx1,:));
nansum(abs(Ts1)>3)
Ts2=BWAS_Tregression([score22(indx2,:),cov_data(indx2,:)],area_data1(indx2,:));
nansum(abs(Ts2)>3)

Ts1(isnan(Ts1))=0;
Ts2(isnan(Ts2))=0;
% 
pv1=normcdf(abs(Ts1),'upper')*2;
pv1(isnan(pv1))=1;
fdrthre1=FDR(pv1,0.05);

pv2=normcdf(abs(Ts2),'upper')*2;
pv2(isnan(pv2))=1;
fdrthre2=FDR(pv2,0.05);

save(['fes_y_ss_fc','_area_withp_indep.mat'],'Ts1','pv1')
save(['pmq_y_ss_mean','_area_withp_indep.mat'],'Ts2','pv2')












