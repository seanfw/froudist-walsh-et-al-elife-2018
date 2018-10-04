% Code to reproduce Figure 6 from Froudist-Walsh et al., eLife, 2018
% (display modules, force-directed graph of network connectivity
% and within-module connectivity changes in modules)
% See readme file for dependencies
% Sean Froudist-Walsh (2018)


close all; clear all; clc;

%% Adjust granularity
lambda = 1;

%% load data
load conn_matrix_no_hipp.mat
load atlas_areas_no_hipp.mat
%% Find most reliable modules

num_iterations = 100; %10000 % Note set to 1 for a quick look, 10000 to assess the most reliable modules
num_regions = size(conn_matrix_mean,1);
num_timepoints = size(conn_matrix_mean,3);

modules_list_mat = nan(num_regions,num_iterations,num_timepoints);
modules_list_pre_reordered = nan(num_regions,num_iterations);
most_common_module_mat = nan(num_regions,num_timepoints);

for current_timepoint = 1:num_timepoints
    for current_iteration = 1:num_iterations

        [modules_list_mat(:,current_iteration,current_timepoint),~] = community_louvain(conn_matrix_mean(:,:,current_timepoint),lambda,[],'negative_asym');
        % give the modules a consistent order
        num_modules = max(modules_list_mat(:,current_iteration,current_timepoint));
        first_module_appearance = nan(num_modules,1);
        modules_mat = nan(num_regions,num_modules);
        
        for current_module = 1:num_modules
            first_module_appearance(current_module) = min(find(modules_list_mat(:,current_iteration,current_timepoint)==current_module));
            modules_mat(:,current_module) = (modules_list_mat(:,current_iteration,current_timepoint)==current_module);

        end
        
        [~,module_order] = sort(first_module_appearance);
    modules_mat_reordered = modules_mat(:,module_order);
    modules_list_mat(:,current_iteration,current_timepoint) = sum(modules_mat_reordered.*repmat(1:num_modules,num_regions,1),2);
    [modules_list, ~, modules_solution] = unique(modules_list_mat(:,:,current_timepoint)','rows');
    most_common_modules = mode(modules_solution);
    most_common_module_mat(:,current_timepoint) = modules_list(most_common_modules,:)';
    end   
end
 
modules_list_pre_no_hipp = most_common_module_mat(:,1);
modules_list_3month_no_hipp = most_common_module_mat(:,2);
modules_list_1year_no_hipp = most_common_module_mat(:,3);

%% Calculate participation coefficient at each timepoint
[part_coeff_pos_pre, part_coeff_neg_pre] = participation_coef_sign(conn_matrix_mean_pre_no_hipp,modules_list_pre_no_hipp);
[part_coeff_pos_3month, part_coeff_neg_3month] = participation_coef_sign(conn_matrix_mean_3month_no_hipp,modules_list_3month_no_hipp);
[part_coeff_pos_1year, part_coeff_neg_1year] = participation_coef_sign(conn_matrix_mean_1year_no_hipp,modules_list_1year_no_hipp);

 
%% reorder timepoint 1 modules to match colours in figure
TCs_L_index = strmatch('TCs_L', region_names_Shen_sequential_no_hipp, 'exact');
PMCdl_R_index = strmatch('PMCdl_R', region_names_Shen_sequential_no_hipp, 'exact');
CCp_R_index= strmatch('CCp_R', region_names_Shen_sequential_no_hipp, 'exact');
TCv_R_index = strmatch('TCv_R', region_names_Shen_sequential_no_hipp, 'exact');

% reorder pre-lesion modules
% parieto-occipital module
modules_list_pre_no_hipp(modules_list_pre_no_hipp==modules_list_pre_no_hipp(CCp_R_index))=7;
% dorsal frontal module
modules_list_pre_no_hipp(modules_list_pre_no_hipp==modules_list_pre_no_hipp(PMCdl_R_index))=8;
% medial temporal module
modules_list_pre_no_hipp(modules_list_pre_no_hipp==modules_list_pre_no_hipp(TCv_R_index))=6;
% anterior temporal/orbitofrontal module
modules_list_pre_no_hipp(modules_list_pre_no_hipp==modules_list_pre_no_hipp(TCs_L_index))=5;

modules_list_pre_no_hipp = modules_list_pre_no_hipp-4;


%% Create nifti files 

% Step 1: load Regional Map atlas
[RM_atlas,scandims,scanscales,scanbpp,scanendian]= read_avw('RM_onMNI_with_hippo_removed.nii.gz');
%
% Step 2: binarise Regional Map atlas 
non_zero_voxels = find(RM_atlas);
RM_atlas_bin = RM_atlas;
RM_atlas_bin(non_zero_voxels) = 1;

%% Create modularity images
%
cd individual_regions

modularity_map_pre = RM_atlas_bin;
hubs_map = RM_atlas_bin;


 clear region_atlas_number;

for current_region = 1:num_regions
    % isolate value of the current region in the atlas - WARNING - this
    % will de-zero-pad the numbers. Make sure nothing gets out of order.
    region_atlas_number(current_region) = str2double(region_numbers_cell_no_hipp{current_region});
    
    %sprintf('Updating modularity map with region %s',region_names_no_hipp{region_atlas_number(current_region)})
     
     
    % load in the region mask
    current_region_mask = read_avw(sprintf('RM_onMNI_region_%d.nii.gz',region_atlas_number(current_region)));
    
    % binarise the region mask
    non_zero_voxels = find(current_region_mask);
    current_region_mask_bin = current_region_mask;
    current_region_mask_bin(non_zero_voxels) = 1;
   
    % value the region mask according to it's modularity group
    current_region_mask_scaled_pre = modules_list_pre_no_hipp(current_region)*current_region_mask_bin;    
    
    % Fill this in in the atlas
    modularity_map_pre = modularity_map_pre+ current_region_mask_scaled_pre;
   
end

% centre the atlas around zero
modularity_map_pre = modularity_map_pre - RM_atlas_bin;


save_avw(modularity_map_pre, '../modularity_map_pre_no_hipp.nii.gz' ,'f',scanscales);

cd ..


%% Non-Matlab part - create volume to surface mapping

!echo "warping modules map pre from MNI to F99"
!applywarp -i modularity_map_pre_no_hipp.nii.gz  -o modularity_map_pre_no_hipp_F99 -r F99_struct_brain.nii.gz -w MNI_to_F99_warp6.nii.gz --interp=nn

!echo "create modularity surface pre - left hemisphere"
!wb_command -volume-to-surface-mapping modularity_map_pre_no_hipp_F99.nii.gz lh.fiducial.surf.gii modularity_map_pre_no_hipp_F99_to_LH.func.gii -trilinear
 
!echo "create modularity surface pre - right hemisphere"
!wb_command -volume-to-surface-mapping modularity_map_pre_no_hipp_F99.nii.gz rh.fiducial.surf.gii modularity_map_pre_no_hipp_F99_to_RH.func.gii -trilinear
 
%% Visualise the surfaces - interactive


lh_surf_struct = gifti('lh.fiducial.surf.gii');
rh_surf_struct = gifti('rh.fiducial.surf.gii');

lh_surf_modules_pre = gifti('modularity_map_pre_no_hipp_F99_to_LH.func.gii');
rh_surf_modules_pre = gifti('modularity_map_pre_no_hipp_F99_to_RH.func.gii');

figure('units','normalized','outerposition',[0 0.67 0.25 0.33])
set(gcf,'color','w');
mysurf=plot(lh_surf_struct,lh_surf_modules_pre);
direction1 = [0 1 0];
direction2 = [0 0 1];
rotate(mysurf,direction1,75);
rotate(mysurf,direction2,90);
colormap jet
title('pre-lesion modules')
figure('units','normalized','outerposition',[0 0.3 0.25 0.33])
mysurf=plot(lh_surf_struct,lh_surf_modules_pre);
direction1 = [0 1 0];
direction2 = [0 0 1];
rotate(mysurf,direction1,270);
rotate(mysurf,direction2,270);
colormap jet
title('pre-lesion modules, right hemisphere')

%% Created forced directed layout - both hemispheres, no hippocampus

desired_threshold = 0.75; % chosen for visualisation purposes

% find strongest connections
[~,index] = sort(conn_matrix_mean_pre_no_hipp(:));

% threshold conn mat
num_remaining_conns = round(desired_threshold*length(index));

threshold = conn_matrix_mean_pre_no_hipp(index(num_remaining_conns));
remaining = conn_matrix_mean_pre_no_hipp>threshold;

conn_matrix_mean_pre_thresh = conn_matrix_mean_pre_no_hipp.*remaining;


% find strongest connections
[~,index] = sort(conn_matrix_mean_3month_no_hipp(:));

% threshold conn mat
num_remaining_conns = round(desired_threshold*length(index));

threshold = conn_matrix_mean_3month_no_hipp(index(num_remaining_conns));
remaining = conn_matrix_mean_3month_no_hipp>threshold;

conn_matrix_mean_3month_thresh = conn_matrix_mean_3month_no_hipp.*remaining;


% find strongest connections
[~,index] = sort(conn_matrix_mean_1year_no_hipp(:));

% threshold conn mat
num_remaining_conns = round(desired_threshold*length(index));

threshold = conn_matrix_mean_1year_no_hipp(index(num_remaining_conns));
remaining = conn_matrix_mean_1year_no_hipp>threshold;

conn_matrix_mean_1year_thresh = conn_matrix_mean_1year_no_hipp.*remaining;


conn_matrix_mean_pre_thresh = round(conn_matrix_mean_pre_thresh,4);
G_pre = graph(conn_matrix_mean_pre_thresh, region_names_Shen_sequential_no_hipp,'OmitSelfLoops');

conn_matrix_mean_3month_thresh = round(conn_matrix_mean_3month_thresh,4);
G_3month = graph(conn_matrix_mean_3month_thresh, region_names_Shen_sequential_no_hipp,'OmitSelfLoops');

conn_matrix_mean_1year_thresh = round(conn_matrix_mean_1year_thresh,4);
G_1year = graph(conn_matrix_mean_1year_thresh, region_names_Shen_sequential_no_hipp,'OmitSelfLoops');

figure('units','normalized','outerposition',[0.23 0.6 0.66 0.33])
set(gcf,'color','w');

subplot(1,3,1)

h_pre = plot(G_pre,'layout','force','EdgeColor',[0.5,0.5,0.5],'MarkerSize',10,'LineWidth',0.1,'NodeLabel',[]);
%for i=1:length(modules_list_pre_no_hipp)
for i = [25,26]  
text(h_pre.XData(i)+rand*0.1,h_pre.YData(i)+rand*0.1,region_names_Shen_sequential_no_hipp(i),'fontsize',16);
end
title('pre-lesion','FontSize',20)
caxis([0,max(modules_list_pre_no_hipp)])
G_pre.Nodes.NodeColors = modules_list_pre_no_hipp;
h_pre.NodeCData = G_pre.Nodes.NodeColors;


subplot(1,3,2)
h_3month = plot(G_3month,'layout','force','EdgeColor',[0.5,0.5,0.5],'MarkerSize',10,'LineWidth',0.1,'NodeLabel',[]);
%for i=1:length(modules_list_pre_no_hipp)
for i = [25,26]   
text(h_3month.XData(i)+rand*0.1,h_3month.YData(i)+rand*0.1,region_names_Shen_sequential_no_hipp(i),'fontsize',16);
end
title('3 months post-lesion','FontSize',20)
caxis([0,max(modules_list_pre_no_hipp)])
G_3month.Nodes.NodeColors = modules_list_pre_no_hipp;
h_3month.NodeCData = G_3month.Nodes.NodeColors;

subplot(1,3,3)
h_1year = plot(G_1year,'layout','force','EdgeColor',[0.5,0.5,0.5],'MarkerSize',10,'LineWidth',0.1,'NodeLabel',[]);

for i = [25,26] 
   text(h_1year.XData(i)-rand*0.1,h_1year.YData(i)-rand*0.1,region_names_Shen_sequential_no_hipp(i),'fontsize',16);
end
 set(gca,'XDir','reverse');
set(gca,'YDir','reverse');
title('12 months post-lesion','FontSize',20)
caxis([0,max(modules_list_pre_no_hipp)])
G_1year.Nodes.NodeColors = modules_list_pre_no_hipp;
h_1year.NodeCData = G_1year.Nodes.NodeColors;
colormap jet

%%
% Calculate within module connectivity

timepoints = {'pre','3month','1year'};
conn_matrices_no_hipp(:,:,1) = conn_matrix_mean_pre_no_hipp;
conn_matrices_no_hipp(:,:,2) = conn_matrix_mean_3month_no_hipp;
conn_matrices_no_hipp(:,:,3) = conn_matrix_mean_1year_no_hipp;
num_modules = max(modules_list_pre_no_hipp);

acute_within_module_connectivity = [];
chronic_within_module_connectivity = [];
Shen_to_modules = [];
%%
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


%% Cateye plots - modules - uses Sam Schwarzkopf's cateye function (see README.txt and https://neuroneurotic.net/2015/09/09/visualising-group-data/)

% acute
[P_acute_module_connectivity,ANOVATAB_acute_module_connectivity,~] = anova1(acute_within_module_connectivity_Shen_order,modules_list_pre_no_hipp,'off');

figure('units','normalized','outerposition',[0 0 0.33 0.33])
colormap jet
my_colormap = colormap;
set(gcf,'color','w');
hold on;
for current_module = 1:num_modules
cateye(acute_within_module_connectivity_Shen_order(modules_list_pre_no_hipp==current_module),current_module,my_colormap(64*(current_module/4),:))
end

title('acute')
ylim([-0.5 0.2])
ylabel('change in intra-modular connectivity')
set(gca,'XTick',1:4)
set(gca,'XTickLabel',{'OFC/ant. temp.','post. temporal','parieto-occipital','dorsal frontal'})
dim = [0.35 0 .3 .3];
str = sprintf('F = %0.4f, p = %0.9f', ANOVATAB_acute_module_connectivity{2,5}, round(P_acute_module_connectivity,9));
a=annotation('textbox',dim,'String',str,'FitBoxToText','on');
a.FontSize = 14;
hline = refline([0 0]);
hline.Color = [0.3 0.3 0.3];
set(hline,'LineStyle',':') 

% chronic
[P_chronic_module_connectivity,ANOVATAB_chronic_module_connectivity,STATS] = anova1(chronic_within_module_connectivity_Shen_order,modules_list_pre_no_hipp,'off');

figure('units','normalized','outerposition',[0.34 0 0.33 0.33])
colormap jet
my_colormap = colormap;
set(gcf,'color','w');
hold on;
for current_module = 1:num_modules
cateye(chronic_within_module_connectivity_Shen_order(modules_list_pre_no_hipp==current_module),current_module,my_colormap(64*(current_module/4),:))
end
title('chronic')
ylim([-0.5 0.2])
ylabel('change in intra-modular connectivity')
set(gca,'XTick',1:4)
set(gca,'XTickLabel',{'OFC/ant. temp.','post. temporal','parieto-occipital','dorsal frontal'})
dim = [.5 0.6 .3 .3];
str = sprintf('F = %0.4f, p = %0.16f', ANOVATAB_chronic_module_connectivity{2,5}, round(P_chronic_module_connectivity,16));
a=annotation('textbox',dim,'String',str,'FitBoxToText','on');
a.FontSize = 14;
hline = refline([0 0]);
hline.Color = [0.3 0.3 0.3];
set(hline,'LineStyle',':') 

