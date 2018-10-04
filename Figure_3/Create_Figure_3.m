% Code to reproduce Figure 3 from Froudist-Walsh et al., eLife, 2018
% (display acute and chronic participation coefficient changes, and 
% calculate relationship with cell densities, hubness and hippocampal
% connectivity)
% See readme file for dependencies
% Sean Froudist-Walsh (2018)


%% Clean up
clc; clear all; close all;

%% load data
load correlated_predictors.mat
load atlas_areas_no_hipp_amyg.mat
load conn_matrix_no_hipp.mat


%% Calculate modules and participation coefficient at current lambda

lambda = 1;

    
[modules_list_pre_no_hipp,~] = community_louvain(conn_matrix_mean_pre_no_hipp,lambda,[],'negative_asym');
[modules_list_3month_no_hipp,~] = community_louvain(conn_matrix_mean_3month_no_hipp,lambda,[],'negative_asym');
[modules_list_1year_no_hipp,~] = community_louvain(conn_matrix_mean_1year_no_hipp,lambda,[],'negative_asym');

num_pre_modules = max(modules_list_pre_no_hipp);
num_3month_modules = max(modules_list_3month_no_hipp);
num_1year_modules = max(modules_list_1year_no_hipp);

% Calculate participation coefficient
[part_coeff_pos_pre, part_coeff_neg_pre] = participation_coef_sign(conn_matrix_mean_pre_no_hipp,modules_list_pre_no_hipp);
[part_coeff_pos_3month, part_coeff_neg_3month] = participation_coef_sign(conn_matrix_mean_3month_no_hipp,modules_list_3month_no_hipp);
[part_coeff_pos_1year, part_coeff_neg_1year] = participation_coef_sign(conn_matrix_mean_1year_no_hipp,modules_list_1year_no_hipp);


% Calculate change in participation coefficient
part_coeff_change_acute = part_coeff_pos_3month -part_coeff_pos_pre;
part_coeff_change_chronic = part_coeff_pos_1year - part_coeff_pos_3month;
part_coeff_change_overall= part_coeff_pos_1year -part_coeff_pos_pre;

% Identify hubs via both node strength & participation coefficient 

[strength_pos_pre] = strengths_und_sign(conn_matrix_mean_pre_no_hipp);

strength_pos_pre_norm = (strength_pos_pre - mean(strength_pos_pre))./std(strength_pos_pre);
part_coeff_pos_pre_norm = (part_coeff_pos_pre - mean(part_coeff_pos_pre))./std(part_coeff_pos_pre);

hubness_matrix_pre = [strength_pos_pre_norm' part_coeff_pos_pre_norm];

[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(hubness_matrix_pre);

hubness_continuous = SCORE(:,1);

% remove amygdala to match atlas
left_amygdala_index_sequential = strmatch('Amyg_L', region_names_Shen_sequential_no_hipp, 'exact');
right_amygdala_index_sequential = strmatch('Amyg_R', region_names_Shen_sequential_no_hipp, 'exact');

part_coeff_change_acute_no_amyg = part_coeff_change_acute;
part_coeff_change_acute_no_amyg([right_amygdala_index_sequential left_amygdala_index_sequential],:)=[];

part_coeff_change_chronic_no_amyg = part_coeff_change_chronic;
part_coeff_change_chronic_no_amyg([right_amygdala_index_sequential left_amygdala_index_sequential],:)=[];

hubness_continuous_no_amyg = hubness_continuous;
hubness_continuous_no_amyg([right_amygdala_index_sequential left_amygdala_index_sequential],:)=[];

correlated_predictors(:,3) = hubness_continuous_no_amyg;

% Feature normalise correlated predictors
mean_predictors_mat  = repmat(mean(correlated_predictors),num_regions_no_hipp_or_amyg,1);
std_predictors_mat  = repmat(std(correlated_predictors),num_regions_no_hipp_or_amyg,1);
correlated_predictors_norm = (correlated_predictors - mean_predictors_mat)./std_predictors_mat;

%% Run stats - acute changes in participation coefficient
% Do stepwise regression.
[b_acute_stepwise,~,~,inmodel_acute_stepwise,stats_acute_stepwise,~,~] = stepwisefit(correlated_predictors_norm,part_coeff_change_acute_no_amyg);

% Calculate r-squared for stepwise regression
stats_acute_stepwise_fstat = stats_acute_stepwise.fstat;
r_squared_acute_stepwise = 1 - (stats_acute_stepwise.SSresid/stats_acute_stepwise.SStotal)';
stats_acute_stepwise_pval = stats_acute_stepwise.pval;
stats_acute_stepwise_TSAT = stats_acute_stepwise.TSTAT;
stats_acute_stepwise_PVAL = stats_acute_stepwise.PVAL;
% Create predicted points 
part_coeff_acute_stepwise_prediction = stats_acute_stepwise.intercept + sum(repmat(inmodel_acute_stepwise.*stats_acute_stepwise.B',num_regions_no_hipp_or_amyg,1).*correlated_predictors_norm,2);

% Decorrelate chronic from acute changes (fit GLM and find residual
% variance in chronic changes not explained by acute changes).
glm_acute_chronic = fitglm(part_coeff_change_acute_no_amyg,part_coeff_change_chronic_no_amyg);
part_coeff_change_chronic_residual = glm_acute_chronic.Residuals.LinearPredictor;

% Compare to random data
num_sims = 1000;
rand_t_vals = nan(num_sims,1);
for current_sim = 1:num_sims
    rand_data    = rand(num_regions_no_hipp_or_amyg,3);
    rand_acute   = rand_data(:,2) - rand_data(:,1);
    rand_chronic = rand_data(:,3) - rand_data(:,2);
    glm_acute_chronic_rand = fitglm(rand_acute,rand_chronic);
    rand_t_vals(current_sim) = table2array(glm_acute_chronic_rand.Coefficients(2,3));

end

 acute_chronic_corrected_p = sum(rand_t_vals<table2array(glm_acute_chronic.Coefficients(2,3)))/num_sims;
%% Run stats - chronic changes in participation coefficient

% Do stepwise regression. 
[b_chronic_stepwise,SE,PVAL,inmodel_chronic_stepwise,stats_chronic_stepwise,NEXTSTEP,HISTORY] = stepwisefit(correlated_predictors_norm,part_coeff_change_chronic_residual);

% Calculate r-squared for stepwise regression
stats_chronic_stepwise_fstat = stats_chronic_stepwise.fstat;
r_squared_chronic_stepwise = 1 - (stats_chronic_stepwise.SSresid/stats_chronic_stepwise.SStotal);
stats_chronic_stepwise_pval = stats_chronic_stepwise.pval;
stats_chronic_stepwise_TSAT = stats_chronic_stepwise.TSTAT;
stats_chronic_stepwise_PVAL = stats_chronic_stepwise.PVAL;

% Create predicted points 
part_coeff_chronic_stepwise_prediction = stats_chronic_stepwise.intercept + sum(repmat(inmodel_chronic_stepwise.*stats_chronic_stepwise.B',num_regions_no_hipp_or_amyg,1).*correlated_predictors_norm,2);


%% Create nifti files 

% Step 1: load Regional Map atlas
[RM_atlas,scandims,scanscales,scanbpp,scanendian]= read_avw('RM_onMNI_with_hippo_and_amyg_removed.nii.gz');
%
% Step 2: binarise Regional Map atlas 
non_zero_voxels = find(RM_atlas);
RM_atlas_bin = RM_atlas;
RM_atlas_bin(non_zero_voxels) = 1;

cd individual_regions

part_coeff_acute_map = RM_atlas_bin;
part_coeff_chronic_map = RM_atlas_bin;

for current_region = 1:num_regions_no_hipp_or_amyg
    % isolate value of the current region in the atlas - WARNING - this
    % will de-zero-pad the numbers. Make sure nothing gets out of order.
    region_atlas_number(current_region) = str2double(region_numbers_cell_no_hipp_or_amyg{current_region});
    
    sprintf('Updating participation coefficient maps with region %s',region_names_Shen_sequential_no_hipp_amyg{current_region})

     
    % load in the region mask
    current_region_mask = read_avw(sprintf('RM_onMNI_region_%d.nii.gz',region_atlas_number(current_region)));
    
    % binarise the region mask
    non_zero_voxels = find(current_region_mask);
    current_region_mask_bin = current_region_mask;
    current_region_mask_bin(non_zero_voxels) = 1;
   
    % value the region mask according to it's neuron or glia density
    current_region_mask_scaled_partcoeff_acute    = part_coeff_change_acute_no_amyg(current_region)*current_region_mask_bin;
    current_region_mask_scaled_partcoeff_chronic  = part_coeff_change_chronic_no_amyg(current_region)*current_region_mask_bin;

    % Fill this in in the atlas
    part_coeff_acute_map   = part_coeff_acute_map  + current_region_mask_scaled_partcoeff_acute;
    part_coeff_chronic_map   = part_coeff_chronic_map  + current_region_mask_scaled_partcoeff_chronic;

end
% centre the atlas around zero
part_coeff_acute_map = part_coeff_acute_map  - RM_atlas_bin;
part_coeff_chronic_map = part_coeff_chronic_map  - RM_atlas_bin;

save_avw(part_coeff_acute_map, '../part_coeff_acute_map.nii.gz' ,'f',scanscales);
save_avw(part_coeff_chronic_map, '../part_coeff_chronic_map.nii.gz' ,'f',scanscales);


cd ..

%% Non-Matlab part - create volume to surface mapping

!echo "warping part_coeff_acute_map from MNI to F99"
!applywarp -i part_coeff_acute_map.nii.gz  -o part_coeff_acute_map_F99 -r F99_struct_brain.nii.gz -w MNI_to_F99_warp6.nii.gz --interp=nn

!echo "warping part_coeff_chronic_map from MNI to F99"
!applywarp -i part_coeff_chronic_map.nii.gz  -o part_coeff_chronic_map_F99 -r F99_struct_brain.nii.gz -w MNI_to_F99_warp6.nii.gz --interp=nn
 


!echo "create part_coeff_acute_map surface - left hemisphere"
!wb_command -volume-to-surface-mapping part_coeff_acute_map_F99.nii.gz lh.fiducial.surf.gii part_coeff_acute_map_F99_to_LH.func.gii -trilinear
  
!echo "create part_coeff_acute_map surface - right hemisphere"
!wb_command -volume-to-surface-mapping part_coeff_acute_map_F99.nii.gz rh.fiducial.surf.gii part_coeff_acute_map_F99_to_RH.func.gii -trilinear

!echo "create part_coeff_chronic_map surface - left hemisphere"
!wb_command -volume-to-surface-mapping part_coeff_chronic_map_F99.nii.gz lh.fiducial.surf.gii part_coeff_chronic_map_F99_to_LH.func.gii -trilinear
  
!echo "create part_coeff_chronic_map surface - right hemisphere"
!wb_command -volume-to-surface-mapping part_coeff_chronic_map_F99.nii.gz rh.fiducial.surf.gii part_coeff_chronic_map_F99_to_RH.func.gii -trilinear


%% Visualise the surfaces - interactive

mkdir figures

lh_surf_struct = gifti('lh.fiducial.surf.gii');

lh_surf_part_coeff_acute     = gifti('part_coeff_acute_map_F99_to_LH.func.gii');
lh_surf_part_coeff_chronic   = gifti('part_coeff_chronic_map_F99_to_LH.func.gii');


myfig = figure('units','normalized','outerposition',[0 0 0.5 0.5]);
axis off;
colormap(b2r(-0.5,0.5))
mysurf=plot(lh_surf_struct,lh_surf_part_coeff_acute);
set(gcf,'color','w');
c=colorbar;
caxis([-0.5 0.5]);
ylabel(c,'acute network participation coefficient change','FontSize',20);
direction1 = [0 1 0];
direction2 = [0 0 1];
rotate(mysurf,direction1,75);
rotate(mysurf,direction2,90);
saveas(myfig,'figures/part_coeff_acute_lateral.eps','epsc');
saveas(myfig,'figures/part_coeff_acute_lateral.fig');
rotate(mysurf,direction1,180);
saveas(myfig,'figures/part_coeff_acute_medial.eps','epsc');
saveas(myfig,'figures/part_coeff_acute_medial.fig');


myfig = figure('units','normalized','outerposition',[0 0 0.5 0.5]);
axis off;
colormap(b2r(-0.5,0.5))
mysurf=plot(lh_surf_struct,lh_surf_part_coeff_chronic);
set(gcf,'color','w');
c=colorbar;
caxis([-0.5 0.5]);
ylabel(c,'chronic network participation coefficient change','FontSize',20);
direction1 = [0 1 0];
direction2 = [0 0 1];
rotate(mysurf,direction1,75);
rotate(mysurf,direction2,90);
saveas(myfig,'figures/part_coeff_chronic_lateral.eps','epsc');
saveas(myfig,'figures/part_coeff_chronic_lateral.fig');
rotate(mysurf,direction1,180);
saveas(myfig,'figures/part_coeff_chronic_medial.eps','epsc');
saveas(myfig,'figures/part_coeff_chronic_medial.fig');

%% Plot stats

% Plot scatter with colours

% Get in hipp_pre colours
hipp_pre = correlated_predictors(:,4);
b2r_hipp_scale = b2r(-0.5,0.5);
hipp_conn_color = round(250*(min(hipp_pre + 0.5,1)) + 1);
colors_hipp_pre = b2r_hipp_scale(hipp_conn_color,:);

myfig = figure('units','normalized','outerposition',[0 0 1 0.5]);
set(gcf,'color','w');
subplot(1,2,1)
scatter(part_coeff_change_acute_no_amyg,part_coeff_acute_stepwise_prediction,45,colors_hipp_pre,'filled');
%scatter(part_coeff_change_acute_no_amyg,part_coeff_acute_stepwise_prediction,45,'k','filled');

hold on;
xlabel({'acute change in','participation coefficient'},'FontSize',20);
ylabel({'predicted acute change in','participation coefficient'},'FontSize',20);
h = lsline;
set(h(1),'LineWidth',3,'Color','k');
l = plot([0 0],[-0.7 0.6], '--');
set(l(1),'LineWidth',2,'Color',[0.7,0.7,0.7]);

regions_to_plot = {'PHC_R','TCv_L','PFCdl_R','PFCdm_L'} ;

for current_region = 1:length(regions_to_plot)
current_region_number = strmatch(regions_to_plot{current_region}, region_names_Shen_sequential_no_hipp_amyg, 'exact');
text(part_coeff_change_acute_no_amyg(current_region_number)+0.01,part_coeff_acute_stepwise_prediction(current_region_number),region_names_Shen_sequential_no_hipp_amyg(current_region_number),'fontsize',24);
end
dim = [.14 .6 .3 .3];
str = sprintf('r^2 = %0.2f, p = %0.g ', r_squared_acute_stepwise,stats_acute_stepwise.pval);
a=annotation('textbox',dim,'String',str,'FitBoxToText','on');
a.FontSize = 20;
title('Acute','FontSize',14);
set(gca, 'FontSize', 25);


subplot(1,2,2)
scatter(part_coeff_change_chronic_residual,part_coeff_chronic_stepwise_prediction,45,colors_hipp_pre,'filled')
hold on;
xlabel({'residual chronic change in','participation coefficient'},'FontSize',20);
ylabel({'predicted chronic change in','participation coefficient'},'FontSize',20);
h = lsline;
set(h(1),'LineWidth',3,'Color','k');
l = plot([0 0],[-0.3 0.3], '--');
set(l(1),'LineWidth',2,'Color',[0.7,0.7,0.7]);

regions_to_plot = {'PHC_R','TCv_L','PFCdl_R','PFCdm_L'} ;

for current_region = 1:length(regions_to_plot)
current_region_number = strmatch(regions_to_plot{current_region}, region_names_Shen_sequential_no_hipp_amyg, 'exact');
text(part_coeff_change_chronic_residual(current_region_number)+0.01,part_coeff_chronic_stepwise_prediction(current_region_number),region_names_Shen_sequential_no_hipp_amyg(current_region_number),'fontsize',24);
end

dim = [.58 .6 .3 .3];
str = sprintf('r^2 = %0.2f, p = %0.8f', r_squared_chronic_stepwise, round(stats_chronic_stepwise.pval,8));
a=annotation('textbox',dim,'String',str,'FitBoxToText','on');
a.FontSize = 20;
title('Chronic','FontSize',14);
set(gca, 'FontSize', 25);


saveas(myfig,'figures/part_coeff_glm_scatter_colours.eps','epsc');
saveas(myfig,'figures/part_coeff_glm_scatter_colours.fig');


myfig = figure('units','normalized','outerposition',[0 0.5 0.5 0.5]);
set(gcf,'color','w');
bar(stats_acute_stepwise.TSTAT,'FaceColor',[0.7,0.7,0.7],'EdgeColor',[0,0,0]);
set(gca,'xticklabel',{'neurons','glia etc', 'hubs', 'hipp conn'},'FontSize',20);
ylabel('t statistic');
ylim([-7 5]);
get(gca, 'YTick');
set(gca, 'FontSize', 20);
saveas(myfig,'figures/part_coeff_glm_acute_predictors.eps','epsc');
saveas(myfig,'figures/part_coeff_glm_acute_predictors.fig');

myfig = figure('units','normalized','outerposition',[0.5 0.5 0.5 0.5]);
set(gcf,'color','w');
bar(stats_chronic_stepwise.TSTAT,'FaceColor',[0.7,0.7,0.7],'EdgeColor',[0,0,0]);
set(gca,'xticklabel',{'neurons','glia etc', 'hubs', 'hipp conn'},'FontSize',20);
ylabel('t statistic');
ylim([-7 5]);
get(gca, 'YTick');
set(gca, 'FontSize', 20)
saveas(myfig,'figures/part_coeff_glm_chronic_predictors.eps','epsc');
saveas(myfig,'figures/part_coeff_glm_chronic_predictors.fig');

%% Save modelfit data for comparison with 2 monkey data
part_coeff_change_acute_no_amyg_full_data = part_coeff_change_acute_no_amyg;
part_coeff_change_chronic_no_amyg_full_data = part_coeff_change_chronic_no_amyg;
betas_acute_fulldata = stats_acute_stepwise.B;
betas_chronic_fulldata = stats_chronic_stepwise.B;

save modelfits_fulldata.mat part_coeff_change_acute_no_amyg_full_data part_coeff_change_chronic_no_amyg_full_data betas_acute_fulldata betas_chronic_fulldata correlated_predictors_norm
