% Code to reproduce Figure 2 from Froudist-Walsh et al., eLife, 2018
% (display cell densities, hubness and hippocampal connectivity)
% See readme file for dependencies
% Sean Froudist-Walsh (2018)

%% Load data and atlas

load predictors.mat
load atlas_areas_no_hipp_amyg.mat

%% Create nifti files 

% Step 1: load Regional Map atlas
[RM_atlas,scandims,scanscales,scanbpp,scanendian]= read_avw('RM_onMNI_with_hippo_and_amyg_removed.nii.gz');
%
% Step 2: binarise Regional Map atlas 
non_zero_voxels = find(RM_atlas);
RM_atlas_bin = RM_atlas;
RM_atlas_bin(non_zero_voxels) = 1;


%%
cd individual_regions

neuron_density = predictors(:,1);
glia_density = predictors(:,2);
hubness = predictors(:,3);
hipp_pre = predictors(:,4);

neuron_density_map = RM_atlas_bin;
glia_density_map = RM_atlas_bin;
hubness_map = RM_atlas_bin;
hipp_pre_map = RM_atlas_bin;



%%

for current_region = 1:num_regions_no_hipp_or_amyg
    % isolate value of the current region in the atlas - WARNING - this
    % will de-zero-pad the numbers. Make sure nothing gets out of order.
    region_atlas_number(current_region) = str2double(region_numbers_cell_no_hipp_or_amyg{current_region});
    
    sprintf('Updating neuron and glia density maps with region %s',region_names_Shen_sequential_no_hipp_amyg{current_region})

     
    % load in the region mask
    current_region_mask = read_avw(sprintf('RM_onMNI_region_%d.nii.gz',region_atlas_number(current_region)));
    
    % binarise the region mask
    non_zero_voxels = find(current_region_mask);
    current_region_mask_bin = current_region_mask;
    current_region_mask_bin(non_zero_voxels) = 1;
   
    % value the region mask according to it's neuron or glia density
    current_region_mask_scaled_neuron    = neuron_density(current_region)*current_region_mask_bin;
    current_region_mask_scaled_glia      = glia_density(current_region)*current_region_mask_bin;
    current_region_mask_scaled_hubness   = hubness(current_region)*current_region_mask_bin;
    current_region_mask_scaled_hipp_pre  = hipp_pre(current_region)*current_region_mask_bin;
    
    % Fill this in in the atlas
    neuron_density_map  = neuron_density_map + current_region_mask_scaled_neuron;
    glia_density_map    = glia_density_map + current_region_mask_scaled_glia;
    hubness_map         = hubness_map + current_region_mask_scaled_hubness;
    hipp_pre_map        = hipp_pre_map + current_region_mask_scaled_hipp_pre;
end

% centre the atlas around zero
neuron_density_map = neuron_density_map  - RM_atlas_bin;
glia_density_map = glia_density_map  - RM_atlas_bin;
hubness_map = hubness_map  - RM_atlas_bin;
hipp_pre_map = hipp_pre_map  - RM_atlas_bin;

save_avw(neuron_density_map, '../neuron_density_map.nii.gz' ,'f',scanscales);
save_avw(glia_density_map, '../glia_density_map.nii.gz' ,'f',scanscales);
save_avw(hubness_map, '../hubness_map.nii.gz' ,'f',scanscales);
save_avw(hipp_pre_map, '../hipp_pre_map.nii.gz' ,'f',scanscales);

cd ..

%% Non-Matlab part - create volume to surface mapping

!echo "warping neuron density map from MNI to F99"
!applywarp -i neuron_density_map.nii.gz  -o neuron_density_map_F99 -r F99_struct_brain.nii.gz -w MNI_to_F99_warp6.nii.gz --interp=nn

!echo "warping glia density map from MNI to F99"
!applywarp -i glia_density_map.nii.gz  -o glia_density_map_F99 -r F99_struct_brain.nii.gz -w MNI_to_F99_warp6.nii.gz --interp=nn
 
!echo "warping hubness map from MNI to F99"
!applywarp -i hubness_map.nii.gz  -o hubness_map_F99 -r F99_struct_brain.nii.gz -w MNI_to_F99_warp6.nii.gz --interp=nn

!echo "warping pre-lesion hippocampal connectivity map from MNI to F99"
!applywarp -i hipp_pre_map.nii.gz -o hipp_pre_map_F99 -r F99_struct_brain.nii.gz -w MNI_to_F99_warp6.nii.gz --interp=nn
 
  


!echo "create neuron density surface - left hemisphere"
!wb_command -volume-to-surface-mapping neuron_density_map_F99.nii.gz lh.fiducial.surf.gii neuron_density_map_F99_to_LH.func.gii -trilinear
  
!echo "create neuron density surface - right hemisphere"
!wb_command -volume-to-surface-mapping neuron_density_map_F99.nii.gz rh.fiducial.surf.gii neuron_density_map_F99_to_RH.func.gii -trilinear

!echo "create glia density surface - left hemisphere"
!wb_command -volume-to-surface-mapping glia_density_map_F99.nii.gz lh.fiducial.surf.gii glia_density_map_F99_to_LH.func.gii -trilinear
  
!echo "create glia density surface - right hemisphere"
!wb_command -volume-to-surface-mapping glia_density_map_F99.nii.gz rh.fiducial.surf.gii glia_density_map_F99_to_RH.func.gii -trilinear

!echo "create hubness surface - left hemisphere"
!wb_command -volume-to-surface-mapping hubness_map_F99.nii.gz lh.fiducial.surf.gii hubness_map_F99_to_LH.func.gii -trilinear
  
!echo "create hubness surface - right hemisphere"
!wb_command -volume-to-surface-mapping hubness_map_F99.nii.gz rh.fiducial.surf.gii hubness_map_F99_to_RH.func.gii -trilinear

!echo "create pre-lesion hippocampal connectivity surface - left hemisphere"
!wb_command -volume-to-surface-mapping hipp_pre_map_F99.nii.gz lh.fiducial.surf.gii hipp_pre_map_F99_to_LH.func.gii -trilinear
  
!echo "create pre-lesion hippocampal connectivity surface - right hemisphere"
!wb_command -volume-to-surface-mapping hipp_pre_map_F99.nii.gz rh.fiducial.surf.gii hipp_pre_map_F99_to_RH.func.gii -trilinear


%% Visualise the surfaces - interactive

mkdir figures

lh_surf_struct = gifti('lh.fiducial.surf.gii');

lh_surf_neuron_density = gifti('neuron_density_map_F99_to_LH.func.gii');
lh_surf_glia_density   = gifti('glia_density_map_F99_to_LH.func.gii');
lh_surf_hubness        = gifti('hubness_map_F99_to_LH.func.gii');
lh_surf_hipp_pre       = gifti('hipp_pre_map_F99_to_LH.func.gii');

whiteredmap = [ones(101,1), linspace(1,0,101)', linspace(1,0,101)'];

myfig = figure('units','normalized','outerposition',[0 0.5 0.5 0.5]);
axis off;
colormap(whiteredmap)
mysurf=plot(lh_surf_struct,lh_surf_neuron_density);
set(gcf,'color','w');
caxis([4.3986e+07, 8.4902e+07]);
c=colorbar;
ylabel(c,'neuronal density','FontSize',20);
direction1 = [0 1 0];
direction2 = [0 0 1];

rotate(mysurf,direction1,75);
rotate(mysurf,direction2,90);
saveas(myfig,'figures/neuronal_density_lateral.eps','epsc');
saveas(myfig,'figures/neuronal_density_lateral.fig');
rotate(mysurf,direction1,180);
saveas(myfig,'figures/neuronal_density_medial.eps','epsc');
saveas(myfig,'figures/neuronal_density_medial.fig');

myfig = figure('units','normalized','outerposition',[0.5 0.5 0.5 0.5]);
axis off;
colormap(whiteredmap)
mysurf=plot(lh_surf_struct,lh_surf_glia_density);
set(gcf,'color','w');
caxis([8.3987e+07,   1.0038e+08]);
c=colorbar;
ylabel(c,'non-neuronal cell density','FontSize',20);
direction1 = [0 1 0];
direction2 = [0 0 1];
rotate(mysurf,direction1,75);
rotate(mysurf,direction2,90);
saveas(myfig,'figures/glial_density_lateral.eps','epsc');
saveas(myfig,'figures/glial_density_lateral.fig');
rotate(mysurf,direction1,180);
saveas(myfig,'figures/glial_density_medial.eps','epsc');
saveas(myfig,'figures/glial_density_medial.fig');

myfig = figure('units','normalized','outerposition',[0 0 0.5 0.5]);
axis off;
colormap(b2r(-2,2))
mysurf=plot(lh_surf_struct,lh_surf_hubness);
set(gcf,'color','w');
c=colorbar;
ylabel(c,'hubness','FontSize',20);
direction1 = [0 1 0];
direction2 = [0 0 1];
rotate(mysurf,direction1,75);
rotate(mysurf,direction2,90);
saveas(myfig,'figures/hubness_lateral.eps','epsc');
saveas(myfig,'figures/hubness_lateral.fig');
rotate(mysurf,direction1,180);
saveas(myfig,'figures/hubness_medial.eps','epsc');
saveas(myfig,'figures/hubness_medial.fig');


myfig = figure('units','normalized','outerposition',[0.5 0 0.5 0.5]);
axis off;
colormap(b2r(-0.5,0.5))
mysurf=plot(lh_surf_struct,lh_surf_hipp_pre);
set(gcf,'color','w');
caxis([-0.5 0.5]);
c=colorbar;
ylabel(c,'pre-lesion hippocampal connectivity','FontSize',20);
direction1 = [0 1 0];
direction2 = [0 0 1];
rotate(mysurf,direction1,75);
rotate(mysurf,direction2,90);
saveas(myfig,'figures/hipp_pre_lateral.eps','epsc');
saveas(myfig,'figures/hipp_pre_lateral.fig');
rotate(mysurf,direction1,180);
saveas(myfig,'figures/hipp_pre_medial.eps','epsc');
saveas(myfig,'figures/hipp_pre_medial.fig');