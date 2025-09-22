% Extract network nodal coordinates and radius from TIF stacks
% ------------------------------------------------------------
% Written by Subramanian Sundaram, 7/2020
% Modified by Navaneeth Krishna Rajeeva Pandian, 10/2023 to include
% tortuosity, and radii & node extraction
%
% This script assumes the TIF stacks consist of one channel
% Ch1 is expected to have a surface marker (such as Lectin) or perfused Dextran
%
% References: The code uses the following MATLAB libraries/resources
% Bioformats - https://docs.openmicroscopy.org/bio-formats/5.3.4/users/matlab/index.html
% DRLSE by C. Li- https://www.mathworks.com/matlabcentral/fileexchange/12711-level-set-for-image-segmentation
% Skeleton3D by P. Kollmannsberger - https://www.mathworks.com/matlabcentral/fileexchange/43400-skeleton3d
% Skel2Graph3D by P. Kollmannsberger - https://www.mathworks.com/matlabcentral/fileexchange/43527-skel2graph-3d
%
% Adapted from original version:
% Transient Support from Fibroblasts is Sufficient to Drive Functional
% Vascularization in Engineered Tissues. Advanced Functional Materials (2020) 2003777

close all; clear all; clc     % Close all figures, clear workspace, clear command window

%% Input Setup
% Set path to file to be analyzed
fldr='/Users/...';  % Folder containing the input image stack
filename = 'Dextran.tif';        % Main TIF file to process
filename1 = "Output%s%s%s";      % Output filename template for CSVs
filenamestl = "Output.stl";      % STL output file name for 3D viewing

% Set path to folder with all required libraries
libpath = '/Users/.../Image-reconstruction/lib';
addpath(fullfile(libpath,'bfmatlab'));        % BioFormats
addpath(fullfile(libpath,'DRLSE_v0'));        % DRLSE segmentation
addpath(fullfile(libpath,'skeletonize'));     % 3D skeletonization
addpath(fullfile(libpath,'skel2graph3d'));    % Convert skeleton to graph

fname = fullfile(fldr,filename);      % Full file path to the TIF stack
reader = bfGetReader(fname);          % Read TIF using BioFormats
omeMeta = reader.getMetadataStore();  % Extract OME-TIFF metadata

% Extract image size in pixels
stackSizeX = omeMeta.getPixelsSizeX(0).getValue(); 
stackSizeY = omeMeta.getPixelsSizeY(0).getValue(); 
stackSizeZ = omeMeta.getPixelsSizeZ(0).getValue(); 

% Get voxel dimensions in µm (X, Y, Z)
voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER);
voxelSizeX = voxelSizeX.doubleValue();
voxelSizeY = omeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROMETER);
voxelSizeY = voxelSizeY.doubleValue();
if stackSizeZ ==1
    voxelSizeZ=1;                     % If 2D image, set Z resolution to 1
else
    voxelSizeZ = omeMeta.getPixelsPhysicalSizeZ(0).value(ome.units.UNITS.MICROMETER);
    voxelSizeZ = voxelSizeZ.doubleValue();
end

%% Read and Resample image stack
info=imfinfo(fname);                 % Read image info from the TIF
num_img = stackSizeZ;               % Number of Z-slices
imwid=stackSizeX; imht=stackSizeY;  % Width and height
imgperch = num_img;                 % Assume 1 channel

% Normalize to isotropic voxel dimensions: resize grid to nx x ny x nz
nx=256;
ny=round((voxelSizeY/voxelSizeX)*stackSizeY/stackSizeX*nx);
nz=round((voxelSizeZ/voxelSizeX)*stackSizeZ/stackSizeX*nx);

% Resample volume to isotropic space using grid interpolation
[y,x,z]=ndgrid(linspace(1,imht,ny),linspace(1,imwid,nx),linspace(1,num_img,nz));
strvol = struct('volume',[]);            % Struct to store volume
imageStack=zeros(imht,imwid,imgperch);   % Original image stack

% Compute average voxel unit after resizing (in µm)
voxelUnit=mean([(voxelSizeX*stackSizeX)/nx, ...
                (voxelSizeY*stackSizeY)/ny, ...
                (voxelSizeZ*stackSizeZ)/nz]);

%% Read and Interpolate Image Stack
num_channels=1; 
for j = 1:num_channels
    for k = 1:imgperch
        % Read image slice-by-slice
        currentImage = imread(fname, (j-1)*imgperch+k, 'Info', info);
        imageStack(:,:,k) = currentImage;
    end
    fprintf('<file> done reading stack...\n')

    % 3D interpolation and median filter
    strvol(j).volume=interp3(medfilt3(imageStack),x,y,z); 
    fprintf(['<file> done interp' num2str(j) '...\n'])
end
%% Smoothing and Normalization
fullStack = strvol(1).volume;
fullStack = imgaussfilt3(fullStack,1,'FilterDomain','spatial');   % Gaussian smoothing

coords = [1 1; nx ny];         % ROI defaults to full image
xcoords = coords(:,1); ycoords = coords(:,2);
x_min = min(xcoords); y_min = min(ycoords);
x_max = max(xcoords); y_max = max(ycoords);

fullStack = fullStack/max(fullStack(:));                        % Normalize intensity
fullStack=fullStack(y_min:y_max,x_min:x_max,:);                % Apply ROI
start_slice = 1;

% Optional: auto-pick first slice with intensity threshold
%start_slice=find(squeeze(max(max(fullStack,[],1),[],2))>0.7,1,'first');
fullStack=fullStack(:,:,start_slice:end);   % Remove pre-slices before start
%% Final Preprocessing and Segmentation
for i=1:2
    fullStack = smooth3(fullStack);  % 3D spatial smoothing
end

[vasc_net]=optimize_threshold(fullStack);  % GUI-assisted segmentation (DRLSE)
vasc_net=imfill(vasc_net,'holes');        % Fill enclosed holes

fprintf('\n')
disp(['<Result> % Volume coverage: ' num2str(sum(vasc_net(:))/numel(vasc_net)*100) '[%]'])
Volume1 = sum(vasc_net(:));               % Total vessel volume

%% Remove Small Isolated Components (Noise Cleanup)
% Identify connected components in the binary segmented network
vasc_net_objects = bwconncomp(vasc_net);

% Measure pixel count for each component
vasc_net_numPixels = cellfun(@numel,vasc_net_objects.PixelIdxList);

% Find indices of small objects (less than 500 pixels)
small_obj_list = find(vasc_net_numPixels < 500);

% Remove small disconnected objects from the network
for i = 1:length(small_obj_list)
    vasc_net(vasc_net_objects.PixelIdxList{small_obj_list(i)}) = 0;
end

%% Skeletonization and Surface Visualization
% Extract 3D skeleton from the cleaned binary volume
vasc_skel = Skeleton3D(vasc_net);

% Plot the vascular surface using isosurface + overlay skeleton
figure();
col=[.7 .7 .8];   % Light gray for surface
hiso = patch(isosurface(vasc_net,0),'FaceColor',col,'EdgeColor','none','FaceAlpha',1);
axis equal; axis off;
lighting phong;
isonormals(vasc_net,hiso);
alpha(0.9);
set(gca,'DataAspectRatio',[1 1 1]);
camlight; hold on;

% Get coordinates of skeleton points for plotting
[w1,l1,h1] = size(vasc_skel);
[x,y,z] = ind2sub([w1,l1,h1], find(vasc_skel(:)));
plot3(y,x,z,'square','Markersize',0.25,'MarkerFaceColor','r','Color','r');

title('Vascular network');
set(gcf,'Color','white');
view(140,80);   % Set 3D view angle

%% Graph Conversion of Skeleton (First Pass)
% Convert the skeleton into a graph: nodes and links
[~, node, link] = Skel2Graph3D(vasc_skel, 0);

% Compute total number of link segments
netlength = sum(cellfun('length', {node.links}));

% Convert graph back to skeleton and re-run to ensure convergence
skel_from_graph = Graph2Skel3D(node, link, w1, l1, h1);
[~, node2, link2] = Skel2Graph3D(skel_from_graph, 0);
netlength_new = sum(cellfun('length', {node2.links}));

% Iterate until the number of links converges
while(netlength_new ~= netlength)
    netlength = netlength_new;
    skel_from_graph = Graph2Skel3D(node2, link2, w1, l1, h1);
    [~, node2, link2] = Skel2Graph3D(skel_from_graph, 0);
    netlength_new = sum(cellfun('length', {node2.links}));
end
%% Visualize Graph-Based Skeleton Structure
% Overlay graph links and nodes over the vessel volume
figure; hold on;

% Define colors for terminal and internal nodes
endline = [0,0,0];             % Gray for end branches
endmarker = 1/255*[200,200,200];
nodeline = [0,0.4470,0.7410];  % Blue for internal branches
nodemarker = [0.9290,0.6940,0.1250];

% Loop through all nodes to draw connected links
for i = 1:length(node2)
    x1 = node2(i).comx;
    y1 = node2(i).comy;
    z1 = node2(i).comz;

    for j = 1:length(node2(i).links)
        % Choose color based on whether this or its connected node is terminal
        if node2(i).ep == 1 || node2(node2(i).conn(j)).ep == 1
            col = endline;
        else
            col = nodeline;
        end

        % Draw link segment between consecutive points
        for k = 1:length(link2(node2(i).links(j)).point) - 1
            [x3,y3,z3] = ind2sub([w1,l1,h1], link2(node2(i).links(j)).point(k));
            [x2,y2,z2] = ind2sub([w1,l1,h1], link2(node2(i).links(j)).point(k+1));
            line([y3 y2],[x3 x2],[z3 z2],'Color',col,'LineWidth',2);
        end
    end

    % Plot node as filled marker
    ncol = nodemarker; if node2(i).ep == 1, ncol = endmarker; end
    plot3(y1,x1,z1,'o','Markersize',9,'MarkerFaceColor',ncol,'Color','k');
end

axis image; axis off; set(gcf,'Color','white');
title('Skeleton - Vascular network');
drawnow; view(-20,45);

%% Calculate Link Lengths (3 Ways)
linklengths = zeros(length(link2),1);     % Cumulative length
linklengths1 = zeros(length(link2),1);    % Point-to-point segment length
linklengths2 = zeros(length(link2),1);    % Direct end-to-end length

for i = 1:length(link2)
    % Get voxel indices of each point on the link
    [lx,ly,lz] = ind2sub([w1,l1,h1], link2(i).point);
    
    % Add the connected nodes at both ends
    lx1 = [node2(link2(i).n1).comx, lx, node2(link2(i).n2).comx];
    ly1 = [node2(link2(i).n1).comy, ly, node2(link2(i).n2).comy];
    lz1 = [node2(link2(i).n1).comz, lz, node2(link2(i).n2).comz];

    lps = [lx1; ly1; lz1]';                    % Assemble into 3D coordinates

    % Method 1: Sum of distances between consecutive points
    lengthmat = diag(pdist2(lps,lps),1);
    linklengths(i) = sum(lengthmat(:));       
    
    % Method 2: Direct distance between first and last point
    DistMat = pdist2(lps,lps);
    linklengths2(i) = DistMat(1,end);
    
    % Method 3: Explicit sum using loop
    for j = 2:length(lx1)
        lengthmat(j) = sqrt((lx1(j)-lx1(j-1))^2 + (ly1(j)-ly1(j-1))^2 + (lz1(j)-lz1(j-1))^2);
    end
    linklengths1(i) = sum(lengthmat(:));
end

%% Visualize skeleton only
m = 0;
figure
for i = 1:length(link2)
    for k = 1:length(link2(i).point)-1
        [x1,y1,z1] = ind2sub([w1,l1,h1], link2(i).point(1));
        [x3,y3,z3] = ind2sub([w1,l1,h1], link2(i).point(k));
        [x2,y2,z2] = ind2sub([w1,l1,h1], link2(i).point(k+1));
        
        % Draw line between points in the link
        line([y3 y2],[x3 x2],[z3 z2],'Color',[0,0,0],'LineWidth',2);
        m = m+1;
    end
    
    % Optionally highlight end of link
    line([y1 y2],[x1 x2],[z1 z2],'Color',[1,0,0],'LineWidth',3);
end

axis image; axis off; set(gcf,'Color','white');
title('Skeleton - Vascular network')
drawnow; view(0,90);

%% Link Length Summary and Preparation for Radius Computation
% Get number of voxel points per link (used for weighted averages later)
k_linklens = zeros(length(link2),1);
for i = 1:length(link2)
    k_linklens(i) = length(link2(i).point) + 1;
end
klen = sum(k_linklens);   % Total number of steps across all links

%% Extract Vessel Radius by Distance to Boundary
vasc_net_perim = bwperim(vasc_net,6);             % Binary perimeter of vessel
vasc_skel_for_search = zeros(size(vasc_skel));
vasc_skel_for_search(unique(horzcat(link2.point))) = 1;  % Mark skeleton points

[D,idx] = bwdist(vasc_skel_for_search);           % Distance map to skeleton
vasc_dist = vasc_net_perim .* D;                  % Distance only at perimeter voxels
vasc_skel_inds = double(idx) .* vasc_net_perim;   % Label perimeter voxels with nearest skel pt index

nonzero_inds = find(vasc_skel_inds);
vasc_skel_inds = vasc_skel_inds(nonzero_inds);
vasc_dist1 = vasc_dist(nonzero_inds);
vasc_net_perim1 = find(vasc_net_perim);

%% Group Boundary Voxels by Skeleton Points
% Group perimeter voxels by nearest skeleton index
[grps, skel_inds] = findgroups(vasc_skel_inds);   
bound_dist = accumarray(grps, vasc_dist1, [], @(v){v});  % Distance values per group
Ellipse_1 = splitapply(@(x){x}, vasc_net_perim1, grps);  % Original indices per group

%% Compute Mean Radius per Link
linkradius = zeros(length(link2),1);       % Mean radius per link
linkradlength = 0;
for i = 1:length(link2)
    linkradlength = linkradlength + length(link2(i).point);
end

linkradius_points = zeros(linkradlength,1);
m4 = 0;
for i = 1:length(link2)
    % Match current link's skeleton points to boundary distance groups
    dist_temp = bound_dist(logical(sum(skel_inds == unique(link2(i).point),2)));
    
    % Store radius values per voxel in the link
    for j = 1:length(dist_temp)
        linkradius_points(j+m4) = mean(dist_temp{j});
    end
    m4 = m4 + size(dist_temp);
    linkradius(i) = mean(vertcat(dist_temp{:}));  % Mean radius for entire link
end

%% Remove Dendrites / Short Segments
% Make a copy of the link structure
link3 = link2;

% Initialize array to mark links for removal
links2remove = zeros(length(link2),1);

% Remove links that are either:
% 1. Shorter than 1 × radius
% 2. Have invalid radius (NaN)
for i = 1:length(link2)
    if linklengths(i) < (1 * linkradius(i)) || isnan(linkradius(i))
        links2remove(i) = i;
    end
end

% Remove marked links from both link list and corresponding radii
links2remove = nonzeros(links2remove);
link3(links2remove) = [];
linkradius3 = linkradius;
linkradius3(links2remove) = [];

%% Compute Length and Tortuosity for Cleaned Links
linklengths3 = zeros(length(link3),1);     % Link lengths (filtered)
directdist = zeros(length(link3),1);       % End-to-end distance
Tortuosity_3 = zeros(length(link3),1);     % Tortuosity = length / direct_dist

for i = 1:length(link3)
    [lx,ly,lz] = ind2sub([w1,l1,h1], link3(i).point);
    lx1 = [node2(link3(i).n1).comx, lx, node2(link3(i).n2).comx];
    ly1 = [node2(link3(i).n1).comy, ly, node2(link3(i).n2).comy];
    lz1 = [node2(link3(i).n1).comz, lz, node2(link3(i).n2).comz];

    % Compute total curve length
    lengthmat = zeros(length(lx1)-1,1);
    for j = 2:length(lx1)
        lengthmat(j) = sqrt((lx1(j)-lx1(j-1))^2 + (ly1(j)-ly1(j-1))^2 + (lz1(j)-lz1(j-1))^2);
    end
    linklengths3(i) = sum(lengthmat);

    % Direct distance between endpoints
    directdist(i) = sqrt((lx1(1)-lx1(end))^2 + ...
                         (ly1(1)-ly1(end))^2 + ...
                         (lz1(1)-lz1(end))^2);

    % Tortuosity = actual length / direct distance
    Tortuosity_3(i) = linklengths3(i) / directdist(i);
end

%% Reassign Radii to Filtered Links (pointwise)
% Compute total number of radius points
linkradlength = 0;
for i = 1:length(link3)
    linkradlength = linkradlength + length(link3(i).point);
end

% Match radius values to new link3 structure
linkradius_points3 = zeros(linkradlength,1);
m = 0; n = 0;
for i = 1:length(link2)
    for j = 1:length(link3)
        if (link2(i).n1 == link3(j).n1) && (link2(i).n2 == link3(j).n2)
            for k = 1:length(link2(i).point)
                m = m + 1;
                linkradius_points3(m) = linkradius_points(n+k);
            end
        end
    end
    n = n + length(link2(i).point);
end

%% Plot Filtered Skeleton (link3)
figure
for i = 1:length(link3)
    for k = 1:length(link3(i).point)-1
        [x3,y3,z3] = ind2sub([w1,l1,h1], link3(i).point(k));
        [x2,y2,z2] = ind2sub([w1,l1,h1], link3(i).point(k+1));
        
        % Draw line segment
        line([y3 y2],[x3 x2],[z3 z2],'Color',[0,0,0],'LineWidth',2);
        
        m = m + 1;
        if i < 2  % optional text for debug
            text(y2 + 0.1, x2 + 0.1, [num2str(m)]);
        end
    end
end

%% Divide Tortuous Links into Segments (for Normal Computation)
link4_len = 0;
Tcheck = 1.15;  % Tortuosity threshold for subdivision

% Estimate length of final link list
for i = 1:length(Tortuosity_3)
    if Tortuosity_3(i) > Tcheck
        link4_len = link4_len + 4;  % Subdivide tortuous links
    end
    link4_len = link4_len + 1;
end

link4_conn = zeros(link4_len,7);      % Store (x1,y1,z1,x2,y2,z2,index)
linkradius4 = zeros(link4_len,1);     % Store radius per link
link4_len = 0; newlens = 0; k = 0;     % Initialize counters
len_cutoff = 12;                      % Small links get fewer subdivisions

for i = 1:length(Tortuosity_3)
    if Tortuosity_3(i) <= Tcheck
        % Save as-is
        link4_len = link4_len + 1;
        [x3,y3,z3] = ind2sub([w1,l1,h1], link3(i).point(1));
        [x2,y2,z2] = ind2sub([w1,l1,h1], link3(i).point(end));
        link4_conn(link4_len,1:6) = [x3 y3 z3 x2 y2 z2];
        k = k + length(link3(i).point);
        link4_conn(link4_len,7) = k;
        linkradius4(link4_len) = linkradius3(i);
    end

    if Tortuosity_3(i) > Tcheck
        if length(link3(i).point) > len_cutoff
            newlens = floor(length(link3(i).point)/5);
            for j = 1:5
                link4_len = link4_len + 1;
                if j < 5
                    [x3,y3,z3] = ind2sub([w1,l1,h1], link3(i).point(1+newlens*(j-1)));
                    [x2,y2,z2] = ind2sub([w1,l1,h1], link3(i).point(newlens*j));
                else
                    [x2,y2,z2] = ind2sub([w1,l1,h1], link3(i).point(end));
                    [x3,y3,z3] = ind2sub([w1,l1,h1], link3(i).point(1+newlens*(j-1)));
                end
                link4_conn(link4_len,1:6) = [x3 y3 z3 x2 y2 z2];
                k = k + (j < 5)*newlens + (j == 5)*(length(link3(i).point)-newlens*(j-1));
                link4_conn(link4_len,7) = k;
                linkradius4(link4_len) = linkradius3(i);
            end
        else
            newlens = floor(length(link3(i).point)/2);
            for j = 1:2
                link4_len = link4_len + 1;
                if j < 2
                    [x3,y3,z3] = ind2sub([w1,l1,h1], link3(i).point(1+newlens*(j-1)));
                    [x2,y2,z2] = ind2sub([w1,l1,h1], link3(i).point(newlens*j));
                else
                    [x2,y2,z2] = ind2sub([w1,l1,h1], link3(i).point(end));
                    [x3,y3,z3] = ind2sub([w1,l1,h1], link3(i).point(1+newlens*(j-1)));
                end
                link4_conn(link4_len,1:6) = [x3 y3 z3 x2 y2 z2];
                k = k + (j < 2)*newlens + (j == 2)*(length(link3(i).point)-newlens*(j-1));
                link4_conn(link4_len,7) = k;
                linkradius4(link4_len) = linkradius3(i);
            end
        end
    end
end

%% Replot Full Skeleton + Overlay Segment Subdivisions
figure
% Replot all original links (link2) in black
for i = 1:length(link2)
    for k = 1:length(link2(i).point) - 1
        [x1,y1,z1] = ind2sub([w1,l1,h1], link2(i).point(1));
        [x3,y3,z3] = ind2sub([w1,l1,h1], link2(i).point(k));
        [x2,y2,z2] = ind2sub([w1,l1,h1], link2(i).point(k+1));
        line([y3 y2],[x3 x2],[z3 z2],'Color',[0,0,0],'LineWidth',2);
        m = m + 1;
    end
end

hold on
% Overlay link4 segments (subdivided links based on tortuosity)
for i = 1:size(link4_conn,1)
    line([link4_conn(i,2) link4_conn(i,5)], ...
         [link4_conn(i,1) link4_conn(i,4)], ...
         [link4_conn(i,3) link4_conn(i,6)], ...
         'Color',[0,1,0],'LineWidth',3); % green lines
end

axis image; axis off; set(gcf,'Color','white');
title('Skeleton - Vascular network')
drawnow; view(0,90);

%% Reorder Coordinates for Compatibility

% Switch x-y and y-x (MATLAB image convention vs Cartesian convention)
v = link4_conn(:,1); link4_conn(:,1) = link4_conn(:,2); link4_conn(:,2) = v;
v = link4_conn(:,4); link4_conn(:,4) = link4_conn(:,5); link4_conn(:,5) = v;

%% Build NodeCon_Radii_Table: Segment-Level Summary Table
NodeCon_Radii_Table = zeros(link4_len, 7);  % Columns: [x1,y1,z1,x2,y2,z2,a,b,angle,radius]
for i = 1:link4_len
    NodeCon_Radii_Table(i,1:6) = link4_conn(i,1:6);          % Start & end point
    NodeCon_Radii_Table(i,7) = linkradius4(i);              % Radius
end

%% Histogram of Link Lengths
linklengths = linklengths * voxelUnit;   % Convert to microns
figure
histogram(linklengths,'BinWidth',20)
title('Histogram of lengths (Vascular network)')
xlabel('Length ({\mu}m)')
ylabel('Count (#)')
set(gcf,'Color','white');

%% Compute Radius Metrics (Mean, Weighted)
link_radius_all = linkradius(~isnan(linkradius)) * voxelUnit;
length_wt = linklengths / sum(linklengths);
link_radius_all1 = linkradius .* length_wt;
link_radius_all1 = link_radius_all1(~isnan(link_radius_all1));
link_radius_all1 = sum(link_radius_all1) * voxelUnit;

length_wt1 = linklengths3 / sum(linklengths3);
link_radius_all2 = linkradius3' * length_wt1 * voxelUnit;

%% Convert Table to Microns (except Orientation)
for i = 1:size(NodeCon_Radii_Table,2)
    NodeCon_Radii_Table(:,i) = NodeCon_Radii_Table(:,i) * voxelUnit;
end

%% Save NodeCon_Radii_Table to CSV File
NodeCon_Radii_Table = round(NodeCon_Radii_Table,2);            % Round to 2 decimals
%NodeCon_Radii_Table_links = round(NodeCon_Radii_Table_links,2);

filename2 = '_';
filename3 = datetime("today");
filename4 = '.csv';
Pointtable_name = sprintf(filename1, filename2, filename3, filename4);
writematrix(NodeCon_Radii_Table(:,:), Pointtable_name);

%% Histogram of Vessel Diameters
figure
histogram(2 * link_radius_all,'BinWidth',10)
title('Histogram of vessel diameters (Vascular network)')
xlabel('Diameter ({\mu}m)')
ylabel('Count (#)')
set(gcf,'Color','white');

%% GUI Viewer: Overlay Vessel on Frame Slider
f_final = figure;
setappdata(f_final,'im1',fullStack);               % Original grayscale
setappdata(f_final,'im1_mask',vasc_net);           % Segmented mask

set(f_final,'Color','white');
image1 = imoverlay(fullStack(:,:,1), bwperim(vasc_net(:,:,1)), 'red');
im_handle = imagesc(image1); axis image; axis off;
setappdata(f_final,'mntg',im_handle);

title('Channel 1 - Overlaid')

slider_display = uicontrol('Parent',f_final,'Style','text','Position',[140,2,100,23], ...
    'String',['Frame No. = ','1'],'BackgroundColor','white');
uicontrol('Parent',f_final,'Style','slider','Position',[80,25,220,23], ...
    'value',1, 'min',1, 'max',size(vasc_net,3), ...
    'callback', {@final_frame_update,slider_display}, ...
    'SliderStep',[1/(num_img-1) 1]);
uicontrol('Parent',f_final,'Style','text','Position',[60,21,23,23], ...
    'String','1','BackgroundColor','white');
uicontrol('Parent',f_final,'Style','text','Position',[300,21,23,23], ...
    'String',num2str(size(vasc_net,3)),'BackgroundColor','white');

%% Display Quantifications in Command Window
v_volume = numel(vasc_net) * voxelUnit^3 * 1e-9; % Convert voxel volume to mm^3

disp(['<Result> Average length of vessels: ' num2str(mean(linklengths)) ' [um]'])
disp(['<Result> Average diameter of vessels: ' num2str(mean(2 * link_radius_all)) ', [um]'])
disp(['<Result> Weighted Average diameter of vessels: ' num2str(2 * sum(link_radius_all1)) ', [um]'])
disp(['<Result> Density of vessels per unit volume: ' num2str(length(link2)/v_volume) ' [#/mm^3]'])

term_nodes = sum(cellfun(@sum, {node2.ep}));
total_nodes = length(node2);
disp(['<Result> Branching nodes: ' num2str(total_nodes - term_nodes) '/' num2str(total_nodes)])
disp(['<Result> Terminal nodes: ' num2str(term_nodes) '/' num2str(total_nodes)])
fprintf('--------------------------------------------\n\n')

%% final_frame_update
% Updates the image and mask when the GUI slider is moved. 
function final_frame_update(es, ~, slider_display)
im1 = getappdata(es.Parent, 'im1');            % Original grayscale stack
im1_mask = getappdata(es.Parent, 'im1_mask');  % Segmented mask
h1 = getappdata(es.Parent, 'mntg');            % Image handle

% Update display text to show current slice
slider_display.String = ['Frame Number = ', num2str(round(es.Value))];
frame_no = round(es.Value);                    % Current Z-slice

% Overlay red outline of vessels over grayscale image
image1 = imoverlay(im1(:,:,frame_no), bwperim(im1_mask(:,:,frame_no)), 'red');
set(h1, 'CData', image1);                      % Update displayed image
end

%% optimize_threshold
% Interactive GUI for selecting a good binary threshold and visualizing DRLSE segmentation.
function [network] = optimize_threshold(im1)
num_img = size(im1,3);                   % Number of slices
thresh_ch1 = 0.4;                        % Initial threshold value
im1_mask = imbinarize(im1, thresh_ch1);  % Binary mask with threshold

% Set up GUI figure and store image stacks
f = figure;
setappdata(f, 'im1_orig', im1);          % Original grayscale
setappdata(f, 'im1_work', im1);          % Working copy (modified by contrast adj.)
setappdata(f, 'im1_mask', im1_mask);     % Current binary mask
setappdata(f, 'thresh1', thresh_ch1);
setappdata(f, 'img_network', zeros(size(im1)));  % Final network mask

% Precompute image enhancements (can switch later)
im1_imadjust = zeros(size(im1));
im1_ahe = zeros(size(im1));
im1_slicehe = zeros(size(im1));

for i=1:size(im1,3)
    im1_imadjust(:,:,i)=imadjust(im1(:,:,i));
    im1_ahe(:,:,i)=adapthisteq(im1(:,:,i));
    im1_slicehe(:,:,i)=imhistmatch(im1(:,:,i),im1(:,:,1));
end

% Save these versions for GUI toggling
setappdata(f, 'im1_imadjust', im1_imadjust);
setappdata(f, 'im1_ahe', im1_ahe);
setappdata(f, 'im1_slicehe', im1_slicehe);
set(f,'Position',[800 160 1000 700]);  % Position the figure
set(f,'Color','white');
curFrame = 1;

axes('Parent',f,'position',[0.13 0.2 0.74 0.75]);
subplot(2,3,1);
h1=imshow(im1(:,:,curFrame));
colormap(gray); axis image; axis off
title('Channel 1')


subplot(2,3,2);
h2=imagesc(im1_mask(:,:,curFrame));
title('Channel 1 - Threshold/Init')

axis image; axis off

subplot(2,3,3);
drlse_mask=drlse_optim(im1,im1_mask,curFrame);
h3=imagesc(imoverlay(im1(:,:,curFrame),bwperim(drlse_mask),'r')); axis off; axis equal; colormap(gray);
axis image; axis off
title('Channel 1 - DRLSE Overlay')

setappdata(f,'h1',h1);
setappdata(f,'h2',h2);
setappdata(f,'h3',h3);

% Frame navigation slider
bgcolor = f.Color;
fr_display=uicontrol('Parent',f,'Style','text','Position',[140,2,100,23],...
    'String',['Frame No. = ',num2str(curFrame)],'BackgroundColor',bgcolor);
setappdata(f,'fr_slider_disp',fr_display);
fr_slider = uicontrol('Parent',f,'Style','slider','Position',[80,25,220,23],...
    'value',curFrame, 'min',1, 'max',size(im1,3),'callback',...
    @frame_update, 'SliderStep',[1/(num_img-1) 1]);
uicontrol('Parent',f,'Style','text','Position',[60,21,23,23],...
    'String','1','BackgroundColor',bgcolor);
uicontrol('Parent',f,'Style','text','Position',[300,21,23,23],...
    'String',num2str(size(im1,3)),'BackgroundColor',bgcolor);

% Threshold slider
thch1_display=uicontrol('Parent',f,'Style','text','Position',[460,2,120,23],...
    'String',['T-Channel1 = ',num2str(thresh_ch1)],'BackgroundColor',bgcolor);
setappdata(f,'thch1_slider_disp',thch1_display);
uicontrol('Parent',f,'Style','slider','Position',[400,25,220,23],...
    'value',thresh_ch1, 'min',0, 'max',1,'callback',...
    {@theta1_update,fr_slider}, 'SliderStep',[1/200 1]);
uicontrol('Parent',f,'Style','text','Position',[380,21,23,23],...
    'String','0','BackgroundColor',bgcolor);
uicontrol('Parent',f,'Style','text','Position',[620,21,23,23],...
    'String','1','BackgroundColor',bgcolor);

% Contrast method toggle
bg=uibuttongroup('Parent',f,'Visible','off','Units','pixels','Position',[120,400,760,23],...
    'SelectionChangedFcn',{@newimageadjust,fr_slider});
uicontrol('Parent',bg,'Style','radiobutton','Position',[10,0,150,23],...
    'value',0, 'String','None','HandleVisibility','off')
uicontrol('Parent',bg,'Style','radiobutton','Position',[90,0,150,23],...
    'value',0, 'String','Im-Adjust slices','HandleVisibility','off')
uicontrol('Parent',bg,'Style','radiobutton','Position',[250,0,150,23],...
    'value',0, 'String','Adapt. Hist. Eq.','HandleVisibility','off')
uicontrol('Parent',bg,'Style','radio','Position',[410,0,150,23],...
    'value',0, 'String','Hist. Eq. Slices','HandleVisibility','off')
bg.Visible='on';

% 'Done' button triggers actual segmentation
buttonDone = uicontrol('Parent',f,'Style','pushbutton','String','Done',...
    'Position', [950 25 40 25]);
buttonDone.Callback = @donePushed;
uiwait(f);   % Pause code until 'Done' is clicked
disp('Done optim...')

network = logical(getappdata(f,'img_network'));  % Return final DRLSE network
close(f);
end

%% theta1_update
% Callback function to update threshold via slider.
function theta1_update(es,~,fr_slider)
slider_display = getappdata(es.Parent,'thch1_slider_disp');
slider_display.String = ['T-Channel1 = ', num2str(es.Value)];

im1 = getappdata(es.Parent,'im1_work');
im1_mask = imbinarize(im1, es.Value);               % New binary mask
setappdata(es.Parent, 'im1_mask', im1_mask);        % Update in GUI
setappdata(es.Parent, 'thresh1', es.Value);         % Save value
frame_update(fr_slider, []);                        % Refresh GUI
end

%% newimageadjust
% Callback for contrast enhancement toggle buttons.
function newimageadjust(es, ev, fr_slider)
switch ev.NewValue.String
    case 'None'
        setappdata(es.Parent,'im1_work',getappdata(es.Parent,'im1_orig'));
        setappdata(es.Parent,'im2_work',getappdata(es.Parent,'im2_orig'));
    case 'Im-Adjust slices'
        setappdata(es.Parent,'im1_work',getappdata(es.Parent,'im1_imadjust'));
        setappdata(es.Parent,'im2_work',getappdata(es.Parent,'im2_imadjust'));
    case 'Adapt. Hist. Eq.'
        setappdata(es.Parent,'im1_work',getappdata(es.Parent,'im1_ahe'));
        setappdata(es.Parent,'im2_work',getappdata(es.Parent,'im2_ahe'));
    case 'Hist. Eq. Slices'
        setappdata(es.Parent,'im1_work',getappdata(es.Parent,'im1_slicehe'));
        setappdata(es.Parent,'im2_work',getappdata(es.Parent,'im2_slicehe'));
end

% Update thresholded mask with new image enhancement
im1_mask=imbinarize(getappdata(es.Parent,'im1_work'),getappdata(es.Parent,'thresh1'));
setappdata(es.Parent,'im1_mask',im1_mask);
im2_mask=imbinarize(getappdata(es.Parent,'im2_work'),getappdata(es.Parent,'thresh2'));
setappdata(es.Parent,'im2_mask',im2_mask);
frame_update(fr_slider,[]);
end

%% drlse_optim
% Core function that applies Distance Regularized Level Set Evolution (DRLSE) to one slice.
function drlse_bw = drlse_optim(Img, Img_bw, frame_no)
Img = Img(:,:,frame_no) * 255;            % Normalize grayscale slice

% DRLSE parameters
timestep = 1;
mu = 0.2 / timestep;
iter_inner = 5;
iter_outer = 20;
lambda = 0;
alfa = 0;
epsilon = 1.5;

% Edge-stopping function based on image gradient
sigma = .8;
G = fspecial('gaussian',15,sigma);
Img_smooth = conv2(Img, G, 'same');
[Ix, Iy] = gradient(Img_smooth);
f = Ix.^2 + Iy.^2;
g = 1 ./ (1 + f);  % Lower values at edges

% Initialize level set function (phi)
c0 = 2;
phi = -Img_bw(:,:,frame_no)*2*c0 + c0;

% Choose potential function (default: double-well)
potentialFunction = 'double-well';

% Perform DRLSE evolution
for n = 1:iter_outer
    phi = drlse_edge(phi, g, lambda, mu, alfa, epsilon, timestep, iter_inner, potentialFunction);
end

drlse_bw = imbinarize(-phi, 0);   % Final segmentation mask
end

%% Core function that applies Distance Regularized Level Set Evolution (DRLSE) to one slice.
function donePushed(es, ~)
uiresume(es.Parent);                   % Resume paused GUI
set(es.Parent,'Visible','off');       % Hide figure

% Apply DRLSE on each slice individually
im1 = getappdata(es.Parent,'im1_work');
im1_mask = getappdata(es.Parent,'im1_mask');
Img_network = getappdata(es.Parent,'img_network');

for i = 1:size(im1,3)
    disp(['<DRLSE> processing slice ' num2str(i) '...']);
    Img_network(:,:,i) = drlse_optim(im1, im1_mask, i);
end

setappdata(es.Parent,'img_network', Img_network);  % Save final mask
end

%% frame_update
% Refresh the three subplots in optimize_threshold when navigating frames.
function frame_update(es,~)
im1 = getappdata(es.Parent,'im1_work');
im2 = getappdata(es.Parent,'im1_mask');

h1 = getappdata(es.Parent,'h1');
h2 = getappdata(es.Parent,'h2');
h3 = getappdata(es.Parent,'h3');

slider_display = getappdata(es.Parent,'fr_slider_disp');
slider_display.String = ['Frame Number = ', num2str(round(es.Value))];
set(h1, 'CData', im1(:,:,round(es.Value)));
set(h2, 'CData', im2(:,:,round(es.Value)));
% Apply DRLSE again on current slice and show overlay
set(h3, 'CData', imoverlay(im1(:,:,round(es.Value)), ...
    bwperim(drlse_optim(im1, im2, round(es.Value))), 'r'));
end
