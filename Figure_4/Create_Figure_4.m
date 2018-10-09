% Code to reproduce Figure 4 from Froudist-Walsh et al., eLife, 2018
% (display acute and chronic within-module connectivity changes, and 
% calculate relationship with cell densities, hubness and hippocampal
% connectivity)
% See readme file for dependencies
% Sean Froudist-Walsh (2018)


%% Clean up
clc; clear all; close all;

%% load data
load predictors.mat
load atlas_areas_no_hipp_amyg.mat
load conn_matrix_no_hipp.mat


%% Calculate modules, participation coefficient and hubness at current lambda

lambda = 1;

[modules_list_pre_no_hipp,~] = community_louvain(conn_matrix_mean_pre_no_hipp,lambda,[],'negative_asym');
[modules_list_3month_no_hipp,~] = community_louvain(conn_matrix_mean_3month_no_hipp,lambda,[],'negative_asym');
[modules_list_1year_no_hipp,~] = community_louvain(conn_matrix_mean_1year_no_hipp,lambda,[],'negative_asym');

num_pre_modules = max(modules_list_pre_no_hipp);
num_3month_modules = max(modules_list_3month_no_hipp);
num_1year_modules = max(modules_list_1year_no_hipp);

% Calculate participation coefficient 
[part_coeff_pos_pre, part_coeff_neg_pre] = participation_coef_sign(conn_matrix_mean_pre_no_hipp,modules_list_pre_no_hipp);

% Calculate within module connectivity
conn_matrices_no_hipp = nan(size(conn_matrix_mean_pre_no_hipp,1),size(conn_matrix_mean_pre_no_hipp,2),3);
timepoints = {'pre','3month','1year'};
conn_matrices_no_hipp(:,:,1) = conn_matrix_mean_pre_no_hipp;
conn_matrices_no_hipp(:,:,2) = conn_matrix_mean_3month_no_hipp;
conn_matrices_no_hipp(:,:,3) = conn_matrix_mean_1year_no_hipp;
num_modules = max(modules_list_pre_no_hipp);

acute_within_module_connectivity = [];
chronic_within_module_connectivity = [];
Shen_to_modules = [];

for current_module = 1:num_modules
    
    num_nodes_in_module = sum(modules_list_pre_no_hipp==current_module);
    
    if num_nodes_in_module > 1
        
    % Isolate the nodes in the current module
    module_conn_mat = conn_matrices_no_hipp(modules_list_pre_no_hipp==current_module,modules_list_pre_no_hipp==current_module,:);
    
    % Remove the self-connections
    module_conn_mat_no_diag_transpose = permute(module_conn_mat,[2,1,3]);
    module_conn_mat_no_diag_transpose(find(repmat(eye(size(module_conn_mat_no_diag_transpose,1)),1,1,3)))=[];
    module_conn_mat_no_diag = permute(reshape(module_conn_mat_no_diag_transpose,size(module_conn_mat,1)-1,size(module_conn_mat,1),3),[2,1,3]);
    % Get the mean within module connection strength for each node
    intramodular_conn_per_node= squeeze(mean(module_conn_mat_no_diag,2));
    
    % Get the mean within module connection strength across all nodes
    intramodular_strength_mean = mean(intramodular_conn_per_node);

    intramodular_conn_change_per_node_acute = intramodular_conn_per_node(:,2)-intramodular_conn_per_node(:,1);
    intramodular_conn_change_per_node_chronic = intramodular_conn_per_node(:,3)-intramodular_conn_per_node(:,2);

    else
    intramodular_conn_change_per_node_acute = 0;
    intramodular_conn_change_per_node_chronic = 0;
 
        
    end
    
    acute_within_module_connectivity = [acute_within_module_connectivity;intramodular_conn_change_per_node_acute];
    chronic_within_module_connectivity = [chronic_within_module_connectivity;intramodular_conn_change_per_node_chronic];

    Shen_to_modules = [Shen_to_modules;find(modules_list_pre_no_hipp==current_module)];
end

[Y, modules_to_Shen]=sort(Shen_to_modules);
acute_within_module_connectivity_Shen_order = acute_within_module_connectivity(modules_to_Shen);
chronic_within_module_connectivity_Shen_order = chronic_within_module_connectivity(modules_to_Shen);

% Identify hubs via both node strength & participation coefficient

[strength_pos_pre] = strengths_und_sign(conn_matrix_mean_pre_no_hipp);

strength_pos_pre_norm = (strength_pos_pre - mean(strength_pos_pre))./std(strength_pos_pre);
part_coeff_pos_pre_norm = (part_coeff_pos_pre - mean(part_coeff_pos_pre))./std(part_coeff_pos_pre);

hubness_matrix_pre = [strength_pos_pre_norm' part_coeff_pos_pre_norm];

[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(hubness_matrix_pre);

hubness_continuous = SCORE(:,1);

% remove amygdala to match with cell density data
left_amygdala_index_sequential = strmatch('Amyg_L', region_names_Shen_sequential_no_hipp, 'exact');
right_amygdala_index_sequential = strmatch('Amyg_R', region_names_Shen_sequential_no_hipp, 'exact');

acute_within_module_connectivity_Shen_order_no_amyg = acute_within_module_connectivity_Shen_order;
acute_within_module_connectivity_Shen_order_no_amyg([right_amygdala_index_sequential left_amygdala_index_sequential])=[];

chronic_within_module_connectivity_Shen_order_no_amyg = chronic_within_module_connectivity_Shen_order;
chronic_within_module_connectivity_Shen_order_no_amyg([right_amygdala_index_sequential left_amygdala_index_sequential])=[];

hubness_continuous_no_amyg = hubness_continuous;
hubness_continuous_no_amyg([right_amygdala_index_sequential left_amygdala_index_sequential],:)=[];

predictors(:,3) = hubness_continuous_no_amyg;

% Feature normalise the data.
mean_predictors_mat  = repmat(mean(predictors),num_regions_no_hipp_or_amyg,1);
std_predictors_mat  = repmat(std(predictors),num_regions_no_hipp_or_amyg,1);
predictors_norm = (predictors - mean_predictors_mat)./std_predictors_mat;

%% Run stats - acute changes in within-module connectivity

% Do stepwise regression.
[b_acute_stepwise,~,~,inmodel_acute_stepwise,stats_acute_stepwise,~,~]   = stepwisefit(predictors_norm,acute_within_module_connectivity_Shen_order_no_amyg);
% Calculate r-squared for stepwise regression
stats_acute_stepwise_fstat = stats_acute_stepwise.fstat;
r_squared_acute_stepwise = 1 - (stats_acute_stepwise.SSresid/stats_acute_stepwise.SStotal)';
stats_acute_stepwise_pval = stats_acute_stepwise.pval;
stats_acute_stepwise_TSAT = stats_acute_stepwise.TSTAT;
stats_acute_stepwise_PVAL = stats_acute_stepwise.PVAL;

% Create predicted points 
within_module_acute_stepwise_prediction = stats_acute_stepwise.intercept + sum(repmat(inmodel_acute_stepwise.*stats_acute_stepwise.B',num_regions_no_hipp_or_amyg,1).*predictors_norm,2);

% Decorrelate chronic from acute changes (fit GLM and find residual
% variance in chronic changes not explained by acute changes).
glm_acute_chronic = fitglm(acute_within_module_connectivity_Shen_order_no_amyg,chronic_within_module_connectivity_Shen_order_no_amyg);
within_module_change_chronic_residual = glm_acute_chronic.Residuals.LinearPredictor;

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

acute_chronic_corrected_p = sum(rand_t_vals>table2array(glm_acute_chronic.Coefficients(2,3)))/num_sims;
%% Run stats - chronic changes in within-module connectivity

% Do stepwise regression. 
[b_chronic_stepwise,SE,PVAL,inmodel_chronic_stepwise,stats_chronic_stepwise,NEXTSTEP,HISTORY] = stepwisefit(predictors_norm,within_module_change_chronic_residual);
stats_chronic_stepwise_fstat = stats_chronic_stepwise.fstat;
r_squared_chronic_stepwise = 1 - (stats_chronic_stepwise.SSresid/stats_chronic_stepwise.SStotal)';
stats_chronic_stepwise_pval = stats_chronic_stepwise.pval;
stats_chronic_stepwise_TSAT = stats_chronic_stepwise.TSTAT;
stats_chronic_stepwise_PVAL = stats_chronic_stepwise.PVAL;

%Create predicted points 
within_module_chronic_stepwise_prediction = stats_chronic_stepwise.intercept + sum(repmat(inmodel_chronic_stepwise.*stats_chronic_stepwise.B',num_regions_no_hipp_or_amyg,1).*predictors_norm,2);

%% Create nifti files 

% Step 1: load Regional Map atlas
[RM_atlas,scandims,scanscales,scanbpp,scanendian]= read_avw('RM_onMNI_with_hippo_and_amyg_removed.nii.gz');
%
% Step 2: binarise Regional Map atlas 
non_zero_voxels = find(RM_atlas);
RM_atlas_bin = RM_atlas;
RM_atlas_bin(non_zero_voxels) = 1;

cd individual_regions

within_module_acute_map = RM_atlas_bin;
within_module_chronic_map = RM_atlas_bin;


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
    current_region_mask_scaled_within_module_acute    = acute_within_module_connectivity_Shen_order_no_amyg(current_region)*current_region_mask_bin;
    current_region_mask_scaled_within_module_chronic      = chronic_within_module_connectivity_Shen_order_no_amyg(current_region)*current_region_mask_bin;

    % Fill this in in the atlas
    within_module_acute_map   = within_module_acute_map  + current_region_mask_scaled_within_module_acute;
    within_module_chronic_map   = within_module_chronic_map  + current_region_mask_scaled_within_module_chronic;

end
% centre the atlas around zero
within_module_acute_map = within_module_acute_map  - RM_atlas_bin;
within_module_chronic_map = within_module_chronic_map  - RM_atlas_bin;

save_avw(within_module_acute_map, '../within_module_acute_map.nii.gz' ,'f',scanscales);
save_avw(within_module_chronic_map, '../within_module_chronic_map.nii.gz' ,'f',scanscales);


cd ..

%% Non-Matlab part - create volume to surface mapping

!echo "warping within_module_acute_map from MNI to F99"
!applywarp -i within_module_acute_map.nii.gz  -o within_module_acute_map_F99 -r F99_struct_brain.nii.gz -w MNI_to_F99_warp6.nii.gz --interp=nn

!echo "warping within_module_chronic_map from MNI to F99"
!applywarp -i within_module_chronic_map.nii.gz  -o within_module_chronic_map_F99 -r F99_struct_brain.nii.gz -w MNI_to_F99_warp6.nii.gz --interp=nn
 


!echo "create within_module_acute_map surface - left hemisphere"
!wb_command -volume-to-surface-mapping within_module_acute_map_F99.nii.gz lh.fiducial.surf.gii within_module_acute_map_F99_to_LH.func.gii -trilinear

!echo "create within_module_acute_map surface - right hemisphere"
!wb_command -volume-to-surface-mapping within_module_acute_map_F99.nii.gz rh.fiducial.surf.gii within_module_acute_map_F99_to_RH.func.gii -trilinear

!echo "create within_module_chronic_map surface - left hemisphere"
!wb_command -volume-to-surface-mapping within_module_chronic_map_F99.nii.gz lh.fiducial.surf.gii within_module_chronic_map_F99_to_LH.func.gii -trilinear
  
!echo "create within_module_chronic_map surface - right hemisphere"
!wb_command -volume-to-surface-mapping within_module_chronic_map_F99.nii.gz rh.fiducial.surf.gii within_module_chronic_map_F99_to_RH.func.gii -trilinear




%% Visualise the surfaces - interactive

mkdir figures

lh_surf_struct = gifti('lh.fiducial.surf.gii');

lh_surf_within_module_acute     = gifti('within_module_acute_map_F99_to_LH.func.gii');
lh_surf_within_module_chronic   = gifti('within_module_chronic_map_F99_to_LH.func.gii');


myfig = figure('units','normalized','outerposition',[0 0 0.5 0.5]);
axis off;
colormap(b2r(-0.5,0.5))
mysurf=plot(lh_surf_struct,lh_surf_within_module_acute);
set(gcf,'color','w');
c=colorbar;
caxis([-0.5 0.5]);
ylabel(c,'acute within-module connectivity change','FontSize',20);
direction1 = [0 1 0];
direction2 = [0 0 1];
rotate(mysurf,direction1,75);
rotate(mysurf,direction2,90);
saveas(myfig,'figures/within_module_acute_lateral.eps','epsc');
saveas(myfig,'figures/within_module_acute_lateral.fig');
rotate(mysurf,direction1,180);
saveas(myfig,'figures/within_module_acute_medial.eps','epsc');
saveas(myfig,'figures/within_module_acute_medial.fig');


myfig = figure('units','normalized','outerposition',[0 0 0.5 0.5]);
axis off;
colormap(b2r(-0.5,0.5))
mysurf=plot(lh_surf_struct,lh_surf_within_module_chronic);
set(gcf,'color','w');
c=colorbar;
caxis([-0.5 0.5]);
ylabel(c,'chronic within-module connectivity change','FontSize',20);
direction1 = [0 1 0];
direction2 = [0 0 1];
rotate(mysurf,direction1,75);
rotate(mysurf,direction2,90);
saveas(myfig,'figures/within_module_chronic_lateral.eps','epsc');
saveas(myfig,'figures/within_module_chronic_lateral.fig');
rotate(mysurf,direction1,180);
saveas(myfig,'figures/within_module_chronic_medial.eps','epsc');
saveas(myfig,'figures/within_module_chronic_medial.fig');

%% Plot stats
% Scatter plots coloured by hippo conn
% Get in hipp_pre colours
hipp_pre = predictors(:,4);
b2r_hipp_scale = b2r(-0.5,0.5);
hipp_conn_color = round(250*(min(hipp_pre + 0.5,1)) + 1);
colors_hipp_pre = b2r_hipp_scale(hipp_conn_color,:);

myfig = figure('units','normalized','outerposition',[0 0 1 0.5]);
set(gcf,'color','w');
subplot(1,2,1)
scatter(acute_within_module_connectivity_Shen_order_no_amyg,within_module_acute_stepwise_prediction,45,colors_hipp_pre,'filled');

hold on;
xlabel({'acute change in','within module connectivity'},'FontSize',20);
ylabel({'predicted acute change in','within module connectivity'},'FontSize',20);
h = lsline;
set(h(1),'LineWidth',3,'Color','k');
l = plot([0 0],[-0.1 0.15], '--');
set(l(1),'LineWidth',2,'Color',[0.7,0.7,0.7]);

regions_to_plot = {'PHC_R','TCv_L','PFCdl_R','PFCdm_L'} ;

for current_region = 1:length(regions_to_plot)
current_region_number = strmatch(regions_to_plot{current_region}, region_names_Shen_sequential_no_hipp_amyg, 'exact');
text(acute_within_module_connectivity_Shen_order_no_amyg(current_region_number)+0.01,within_module_acute_stepwise_prediction(current_region_number),region_names_Shen_sequential_no_hipp_amyg(current_region_number),'fontsize',24);
end

dim = [.14 .6 .3 .3];
str = sprintf('r^2 = %0.2f, p = %0.g ', r_squared_acute_stepwise,stats_acute_stepwise.pval);
a=annotation('textbox',dim,'String',str,'FitBoxToText','on');
a.FontSize = 20;
title('Acute','FontSize',14);
set(gca, 'FontSize', 25);



subplot(1,2,2);
scatter(within_module_change_chronic_residual,within_module_chronic_stepwise_prediction,45,colors_hipp_pre,'filled')

hold on;
xlabel({'residual chronic change in','within module connectivity'},'FontSize',20);
ylabel({'predicted chronic change in','within module connectivity'},'FontSize',20);
h = lsline;
set(h(1),'LineWidth',3,'Color','k');
l = plot([0 0],[-0.2 0.2], '--');
set(l(1),'LineWidth',2,'Color',[0.7,0.7,0.7]);

regions_to_plot = {'PHC_R','TCv_L','PFCdl_R','PFCdm_L'} ;

for current_region = 1:length(regions_to_plot)
current_region_number = strmatch(regions_to_plot{current_region}, region_names_Shen_sequential_no_hipp_amyg, 'exact');
text(within_module_change_chronic_residual(current_region_number)+0.01,within_module_chronic_stepwise_prediction(current_region_number),region_names_Shen_sequential_no_hipp_amyg(current_region_number),'fontsize',24);
end

dim = [.58 .6 .3 .3];
str = sprintf('r^2 = %0.2f, p = %0.8f', r_squared_chronic_stepwise, round(stats_chronic_stepwise.pval,8));
a=annotation('textbox',dim,'String',str,'FitBoxToText','on');
a.FontSize = 20;
title('Chronic','FontSize',14);
set(gca, 'FontSize', 25);


saveas(myfig,'figures/within_module_glm_scatter_colours.eps','epsc');
saveas(myfig,'figures/within_module_glm_scatter_colours.fig');


%%

myfig = figure('units','normalized','outerposition',[0 0.5 0.5 0.5]);
set(gcf,'color','w');
bar(stats_acute_stepwise.TSTAT,'FaceColor',[0.7,0.7,0.7],'EdgeColor',[0,0,0]);
set(gca,'xticklabel',{'neurons','glia etc', 'hubs', 'hipp conn'},'FontSize',20);
ylabel('t statistic');
ylim([-7 5]);
get(gca, 'YTick');
set(gca, 'FontSize', 20);
saveas(myfig,'figures/within_module_glm_acute_predictors.eps','epsc');
saveas(myfig,'figures/within_module_glm_acute_predictors.fig');

myfig = figure('units','normalized','outerposition',[0.5 0.5 0.5 0.5]);
set(gcf,'color','w');
bar(stats_chronic_stepwise.TSTAT,'FaceColor',[0.7,0.7,0.7],'EdgeColor',[0,0,0]);
set(gca,'xticklabel',{'neurons','glia etc', 'hubs', 'hipp conn'},'FontSize',20);
ylabel('t statistic');
ylim([-7 5])
get(gca, 'YTick');
set(gca, 'FontSize', 20);
saveas(myfig,'figures/within_module_glm_chronic_predictors.eps','epsc');
saveas(myfig,'figures/within_module_glm_chronic_predictors.fig');


%% Save modelfit data for comparison with 2-monkey data
within_module_acute_full_data = acute_within_module_connectivity_Shen_order_no_amyg;
within_module_chronic_full_data = chronic_within_module_connectivity_Shen_order_no_amyg;
betas_acute_fulldata = stats_acute_stepwise.B;
betas_chronic_fulldata = stats_chronic_stepwise.B;

save modelfits_fulldata.mat within_module_acute_full_data within_module_chronic_full_data betas_acute_fulldata betas_chronic_fulldata predictors_norm
