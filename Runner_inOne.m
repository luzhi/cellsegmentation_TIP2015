function Runner_inOne( DSOption, path_src_nuclei_images_folder, path_src_cyto_images_folder, ...
                       storageCommonPath, storageInitial, storageDist, ...
                       beta_logistic_set, kappa_set, chi_set, loop, ...
                       OUTPUT_PATH)
%%
%==================================
%           Parameters
%==================================
 
%==================================
% Common parameters for GMM-based 
% clump training/testing
%==================================
thr.minCellSize = 1000;
 
%==================================
% Level set paramters for clump 
% refinement (after GMM)
%==================================
iter_in = 2;
iter_out = 3;
alfaGMM = -2.5;
lambdaGMM = 5;
 
%==================================
% Level set parameters for nuclei 
% refinement (after GMM)
%==================================
iter_in_rawNuclei = 10;
iter_out_rawNuclei = 2;
alfa_rawNuclei = 5;
lambda_rawNuclei = 4;

%==================================
%   Level set for segmentation
%==================================
iter_in_extent = 20;
iter_out_extent = 2;

%==================================
%     Radiating distance map
%==================================
% nLinePartition = 10;
min_PtDistance = 5;
min_TinyFragments_DistMap = 500;


%%
%===================================
%        Read source images
%===================================
if strcmp(DSOption, 'EDF')
%+-------------------------------+
%|    Pap smear - EDF Dataset    |
%+-------------------------------+
    %+--------------------------------+
    %|    Source nuclei/cyto image    |
    %+--------------------------------+
    listing_src_edf_img = dir(path_src_nuclei_images_folder);
    imNum_src_nuclei = 0;
    for i = 1:size(listing_src_edf_img,1)
        if listing_src_edf_img(i,1).isdir == 1
            continue;
        end

        if isempty(strfind(listing_src_edf_img(i,1).name, 'GTMask'))
            imNum_src_nuclei = imNum_src_nuclei + 1;
            imNucleiSet{imNum_src_nuclei,1} = imread(strcat(path_src_nuclei_images_folder,listing_src_edf_img(i,1).name));

            if ~eq(size(imNucleiSet{imNum_src_nuclei,1}), [1024, 1024])
                imNucleiSet{imNum_src_nuclei,1} = imresize(imNucleiSet{imNum_src_nuclei,1}, 0.5);  
            end          
        end
    end
    
    %+----------------------------------------+
    %|          For EDF dataset,              |
    %| nuclei and cytoplasm on the same image |
    %+----------------------------------------+
    imCytoSet = imNucleiSet;
    
    %+------------------------------+
    %| Images Number in the dataset |
    %|      is same as either       |
    %|   imNucleiSet or imCytoSet   |
    %+------------------------------+
    imNum = size(imCytoSet,1);
end

%+---------------------------------+
%|	Synthetic Train Dataset        |
%+---------------------------------+
if strcmp(DSOption, 'SyntheticTrain')
    load(strcat(path_src_cyto_images_folder, 'trainset.mat'), 'trainset');
    imCytoSet = trainset;
    imNucleiSet = trainset;
    
    imNum = size(imCytoSet,1);
    
    clear trainset;
end

%+---------------------------------+
%|	Synthetic Test Dataset         |
%+---------------------------------+
if strcmp(DSOption, 'SyntheticTest')
    load(strcat(path_src_cyto_images_folder, 'testset.mat'), 'testset');
    imCytoSet = testset;
    imNucleiSet = testset;
    
    imNum = size(imCytoSet,1);
    
    clear testset;
end
%% 
%=================================
%    Cytoplasm's (clump) Mask
%   by Superpixel & Convex Hull
%=================================
try
    load(strcat(OUTPUT_PATH, storageCommonPath, 'ConvexhullCytoClumpMaskSet.mat'), 'ConvexhullCytoClumpMaskSet');
    fprintf('Step - 1 Raw Clumps boundaries by Convex Hull & Level Set refinement.........done!\n');
catch
    fprintf('Step - 1 Raw Clumps boundaries by Convex Hull & Level Set refinement...\n');
    
    ConvexhullCytoClumpMaskSet = cell(imNum,1);
    
    %+-----------------------------------+
    %| Superpixel Parameters by Data set |
    %+-----------------------------------+
    if strcmp(DSOption, 'EDF')
        spRATIO = 0.5;
        spKERNELSIZE = 2;
        spMAXDIST = 10;
    end
    
    if strcmp(DSOption, 'SyntheticTest') || strcmp(DSOption, 'SyntheticTrain')
        spRATIO = 0.5;
        spKERNELSIZE = 2;
        spMAXDIST = 10;
    end
    
    %+---------------------------------+
    %|              Timer              |
    %+---------------------------------+
    t_preprocessing = zeros(size(imCytoSet,1),1);
    for i = 1:length(imCytoSet)
        tic;
        im_i = imCytoSet{i,1};
        
        fprintf('\tSuperpixel Convex Hull: Image %d\n', i);
        [ ConvexhullCytoClumpMask ] = SuperpixelConvexHull( im_i, spRATIO, spKERNELSIZE, spMAXDIST );
 
        ConvexhullCytoClumpMaskSet{i,1} = logical(ConvexhullCytoClumpMask);
        
        t_preprocessing(i) = toc; 
    end
    
    %+-------------------------------+
    %| Save section variable to file |
    %+-------------------------------+
    save(strcat(OUTPUT_PATH, storageCommonPath, 'ConvexhullCytoClumpMaskSet.mat'),'ConvexhullCytoClumpMaskSet', 't_preprocessing');
end

%%
%====================================
%      Cytoplasm's (clump) Mask
%           by GMM Learning
%====================================
try
    load(strcat(OUTPUT_PATH, storageCommonPath, 'GMMCytoClumpMaskSet.mat'), ...
        'GMMCytoClumpMaskSet', 'gmm_model_cytoClump', 'gmm_model_cytobackground', 't_GMM_train');
    fprintf('Step - 2 GMM-based Training/Testing for accurate clump boundary.........done!\n');
catch
    fprintf('Step - 2 GMM-based Training/Testing for accurate clump boundary...\n');
   
    %+---------------------------------+
    %|              Timer              |
    %+---------------------------------+
    tic;
    
    %+---------------------------------+
    %| Semi-supervised Learning by GMM |
    %+---------------------------------+ 
    idx_TrainImg = (1:imNum)';
    GMMCytoClumpMaskSet = GenerateForeBackgroundMask( ConvexhullCytoClumpMaskSet );
    
    for iter_GMM = 1:10
        if iter_GMM > 1
            GMMCytoClumpMaskSet = GenerateForeBackgroundMask( GMMCytoClumpMaskSet );
        end
        
        %+-------------------------------+
        %| Train GMM model by ALL images |
        %+-------------------------------+
        [ gmm_model_cytoClump, gmm_model_cytobackground ]  = TrainGMMSceneSegmentation( imCytoSet, GMMCytoClumpMaskSet, ...
                                                                                    idx_TrainImg ); 
        %+------------------------------+
        %| Test GMM model by ALL images |
        %+------------------------------+
        for i = 1:imNum
            im_i = imCytoSet{i,1};  
            %+------------------------------+
            %| GMM Posterior for Each Pixel |
            %+------------------------------+
            [ gmm_post ] = TestGMMSceneSegmentation( gmm_model_cytoClump, gmm_model_cytobackground, imCytoSet, i );
            
            [ GMMCytoClumpMask ] = ComputeConfidenceAsScene(gmm_post{1,1}, size(im_i,1), size(im_i,2));
            GMMCytoClumpMaskSet{i,1} = GMMCytoClumpMask; 
        end        
        fprintf('\tGMM Learning Iteration: %d\n', iter_GMM);
    end
   
    t_GMM_train = toc;  
    
    %+-------------------------------+
    %| Save section variable to file |
    %+-------------------------------+
    save(strcat(OUTPUT_PATH, storageCommonPath, 'GMMCytoClumpMaskSet.mat'), ...
        'GMMCytoClumpMaskSet', 'gmm_model_cytoClump', 'gmm_model_cytobackground', 't_GMM_train');
end

%%
%============================
%      Remove Fragments 
%   in GMM Learning Result
%============================
try
    load(strcat(OUTPUT_PATH, storageCommonPath, 'SceneCytoClumpMaskSet.mat'), ...
        'SceneCytoClumpMaskSet', 't_GMM_LevelSetClean');
    fprintf('Step - 3 Clean Noise in GMM Learned Cyto Clump Masks.........done!\n');
catch
    fprintf('Step - 3 Clean Noise in GMM Learned Cyto Clump Masks...\n');
    
    %+---------------------------------+
    %|              Timer              |
    %+---------------------------------+
    t_GMM_LevelSetClean = zeros(imNum,1);
    
    %+--------------------------------------+
    %| Remove Fragments by Level Set Method |
    %+--------------------------------------+    
    SceneCytoClumpMaskSet = cell(imNum,1);
    for i = 1:imNum   
        tic;
        im_i = imCytoSet{i,1};
        GMMCytoClumpMask = GMMCytoClumpMaskSet{i,1};
        
        %+------------------------------------------+
        %| Erase Tiny Fragments by Level Set Method |
        %+------------------------------------------+
        GMMCytoClumpMask = cleanFragmentsLevelSet(im_i, GMMCytoClumpMask, ...
                                                  iter_in, iter_out, alfaGMM, lambdaGMM );
        GMMCytoClumpMask = im2bw(GMMCytoClumpMask);        
        imCBMaskLabels = bwlabel(~GMMCytoClumpMask,8);
        
        %+--------------------------------------------------+
        %| Erase Still Existing Tiny Fragments by Threshold |
        %+--------------------------------------------------+
        for j = 1:max(imCBMaskLabels(:))
            idx = find(imCBMaskLabels == j);
            if size(idx,1) < thr.minCellSize
                GMMCytoClumpMask(idx) = 1;
                fprintf('OK Tiny\n');
            end
        end
        GMMCytoClumpMask = im2bw(GMMCytoClumpMask);
        
        %+-------------------------+
        %| Refined Cyto Clump Mask |
        %+-------------------------+
        SceneCytoClumpMaskSet{i,1} = ~GMMCytoClumpMask;
 
        t_GMM_LevelSetClean(i) = toc;   
    end
    
    %+-------------------------------+
    %| Save section variable to file |
    %+-------------------------------+
    save(strcat(OUTPUT_PATH, storageCommonPath, 'SceneCytoClumpMaskSet.mat'), ...
        'SceneCytoClumpMaskSet', 't_GMM_LevelSetClean');
end
 
%% 
%=================================
% Using MSER to find Nuclei Blobs
%================================= 
try
%     load('xx');
    load(strcat(OUTPUT_PATH, storageCommonPath, 'MSERNucleiBlobSet.mat'), 'MSERNucleiBlobSet');
    disp('Step - 3 To find nuclei.........done!');
catch
    disp('Step - 3 To find nuclei...');
    
    %+-----------------------------+
    %| MSER Parameters by Data set |
    %+-----------------------------+   
    if strcmp(DSOption, 'EDF')
        RegionAreaRange = [120,600];
        ThresholdDelta = 3;
    end

    if strcmp(DSOption, 'SyntheticTest') || strcmp(DSOption, 'SyntheticTrain')
        RegionAreaRange = [120,600];
        ThresholdDelta = 3;
    end
    
    %+----------------+
    %| MSER Blobs Set |
    %+----------------+
    MSERNucleiBlobSet = cell(length(imNucleiSet),1);
 
    %+---------------------------------+
    %|              Timer              |
    %+---------------------------------+
    t_MSER_blobs = zeros(size(imNucleiSet,1),1);
    
    for i = 1:length(imNucleiSet)
        tic;
        I = imNucleiSet{i,1};
 
        % Using MSER to find raw candidate regions for nuclei
        mser_regions = detectMSERFeatures(I, 'RegionAreaRange', RegionAreaRange, 'ThresholdDelta', ThresholdDelta);
        MSERNucleiBlob = false(size(I));
        for j = 1:length(mser_regions)    
            mser_region_xy = mser_regions(j,1).PixelList;
            for k = 1:size(mser_regions(j,1).PixelList,1)
                IIx = mser_region_xy(k,1);
                IIy = mser_region_xy(k,2);
                MSERNucleiBlob(IIy,IIx) = 1;
            end
        end
 
        MSERNucleiBlob = logical(MSERNucleiBlob);
 
        %+-----------------------------+
        %| Denoise by Level Set Method |
        %+-----------------------------+
        MSERNucleiBlob = cleanFragmentsLevelSet(I, MSERNucleiBlob, ... 
                                                iter_in_rawNuclei, iter_out_rawNuclei, ...
                                                alfa_rawNuclei, lambda_rawNuclei );

        MSERNucleiBlob = im2bw(MSERNucleiBlob);
 
        %+----------------+
        %| Form Blobs Set |
        %+----------------+
        MSERNucleiBlobSet{i,1} = MSERNucleiBlob;
 
        t_MSER_blobs(i) = toc;
    end
    
    %+-------------------------------+
    %| Save section variable to file |
    %+-------------------------------+
    save(strcat(OUTPUT_PATH, storageCommonPath, 'MSERNucleiBlobSet.mat'), 'MSERNucleiBlobSet', 't_MSER_blobs');    
end

%% 
%==================================
% Using Rules to find nucleus from 
%       MSERNucleiBlobSet
%================================== 
try
    load(strcat(OUTPUT_PATH, storageCommonPath, 'NucleiMask.mat'), 'NucleiMaskSet');
    fprintf('Step 4 - Filter MSER Nuclei Blobs by Rules...done!\n');
catch
    fprintf('Step 4 - Filter MSER Nuclei Blobs by Rules...');
    NucleiMaskSet = cell(imNum,1);
    
    %+---------------------------------+
    %|      Parameters for Rules       |
    %+---------------------------------+
    if strcmp(DSOption, 'EDF') || strcmp(DSOption, 'SyntheticTest') || strcmp(DSOption, 'SyntheticTrain')
        maxNucleiEccentricity = 0.90;
        minNucleiArea = 100;
        maxNucleiClumpAreaRatio = 0.1;
    end

    %+---------------------------------+
    %|              Timer              |
    %+---------------------------------+
    t_MSER_Clean = zeros(size(imNucleiSet,1),1);
    
    for i = 1:imNum
        tic;
        clumpMask = SceneCytoClumpMaskSet{i,1};
        MSERNucleiBlobMask = ~MSERNucleiBlobSet{i,1};
        nucleiMask_i = zeros(size(MSERNucleiBlobMask));
        
        MSERNucleiBlobStats = regionprops(MSERNucleiBlobMask, 'Eccentricity', 'PixelIdxList', 'Area');
        for j = 1:length(MSERNucleiBlobStats)
            
            %+------------------------------------+
            %| Area Ratio of Nuclei and its Clump |
            %+------------------------------------+
            areaRatio = ComputeAreaRatio_NucleiClump(MSERNucleiBlobStats(j,1).Area, ...
                                                     MSERNucleiBlobStats(j,1).PixelIdxList, clumpMask);
             
            %+-------------------------------------+
            %| Average Intensity of this MSER Blob |
            %+-------------------------------------+
            meanIntensity = mean(imNucleiSet{i,1}(MSERNucleiBlobStats(j,1).PixelIdxList(:)));
            
            %+------------------------------------+
            %|      Filter MSER Nuclei Blobs      |
            %|   by Eccentricity and Area Ratio   |
            %+------------------------------------+
            if MSERNucleiBlobStats(j,1).Area > 300 && meanIntensity >= 200
                continue;
            end
            if MSERNucleiBlobStats(j,1).Eccentricity < maxNucleiEccentricity && ...
               MSERNucleiBlobStats(j,1).Area > minNucleiArea &&...
               areaRatio < maxNucleiClumpAreaRatio && ...
               meanIntensity <= 170
           
                nucleiMask_i( MSERNucleiBlobStats(j,1).PixelIdxList(:) ) = 1;
                
            elseif MSERNucleiBlobStats(j,1).Area > 300 && ...
                   MSERNucleiBlobStats(j,1).Eccentricity < maxNucleiEccentricity && ...
                   areaRatio < maxNucleiClumpAreaRatio && ...
                   meanIntensity <= 200
                
                   nucleiMask_i( MSERNucleiBlobStats(j,1).PixelIdxList(:) ) = 1;
            end
        end
        
        NucleiMaskSet{i,1} = nucleiMask_i;
        
        t_MSER_Clean(i) = toc;
    end

    fprintf('done!\n');
    
    %+-------------------------------+
    %| Save section variable to file |
    %+-------------------------------+
    save(strcat(OUTPUT_PATH, storageCommonPath, 'NucleiMask.mat'), 'NucleiMaskSet', 't_MSER_Clean');
end
 
%% 
%=======================================
% Compute Initial Guess by Distance Map
%=======================================
try
    load('xx');
    fprintf('Step - 5 Compute Initial Guess by Distance Map.........done!\n');
catch
    fprintf('Step - 5 Compute Initial Guess by Distance Map.........\n');
    roughInitialGuessMaskSet_byImage = cell(imNum,1);
    phi_new_masks_set = cell(imNum,1);
    GeoInitialGuessMaskSet_byImage = cell(imNum,1);

    rough_DistMapSetByImage = cell(imNum,1);
    DistMapSetByImage_GeoCentroids = cell(imNum,1);
   
    if strcmp(DSOption, 'EDF')
        imSet = imCytoSet;
    end    

    if strcmp(DSOption, 'SyntheticTest') || strcmp(DSOption, 'SyntheticTrain')
        imSet = imCytoSet;
    end

    for para_id_BETA = 1:size(beta_logistic_set,1)
        beta = beta_logistic_set(para_id_BETA);
        fprintf('beta = %f\n', beta);
        for i = 1:imNum
            %+---------------------------------+
            %|              Timer              |
            %+---------------------------------+
            tic;
            
            fprintf('Image: %d\n', i);
            im_i = imSet{i,1};
            sceneCytoClumpMask_i = SceneCytoClumpMaskSet{i,1};
            nucleiMask_i = logical(NucleiMaskSet{i,1});

            %==================================
            % Get clump & nuclei regions 
            % pixels' index
            %==================================
            sceneCytoClumpStats_i = regionprops(sceneCytoClumpMask_i, 'PixelIdxList');
            nucleiStats_i = regionprops(nucleiMask_i, 'PixelIdxList');
            nucleiCounter = 1;
            
            for j = 1:size(sceneCytoClumpStats_i,1)
                nucleusMask_inside_Clump_j = false(size(im_i));
                clumpMask_j = zeros(size(im_i));
                clumpMask_j( sceneCytoClumpStats_i(j,1).PixelIdxList ) = 1;
                clumpMask_j = logical(clumpMask_j);

                %+------------------------------------+
                %| Find the Nucleus inside this Clump |
                %+------------------------------------+
                insideNucleiNum = 0;
                for n = 1:size(nucleiStats_i,1)                
                    isNucleiInsideClump = isequal(intersect(nucleiStats_i(n,1).PixelIdxList(:),...
                                          sceneCytoClumpStats_i(j,1).PixelIdxList(:)),...
                                          nucleiStats_i(n,1).PixelIdxList(:));
                    if isNucleiInsideClump == 1
                        %+----------------------------+
                        %|   Nuclei Mask of A Clump   |
                        %+----------------------------+
                        nucleusMask_inside_Clump_j( nucleiStats_i(n,1).PixelIdxList ) = 1;
                        insideNucleiNum = insideNucleiNum + 1;
                    end
                end

                %+------------------------------------+
                %| Single Nuclei Means NO Overlapping |
                %| put empty matrix, not involved in  |
                %|       level set evolution          |
                %+------------------------------------+
                if insideNucleiNum == 1
                    nuclei_centroid_xy_i = regionprops(nucleusMask_inside_Clump_j, 'Centroid');
                    roughInitialGuessMaskSet_byImage{i,1}{nucleiCounter,1} = logical(clumpMask_j);
                    roughInitialGuessMaskSet_byImage{i,1}{nucleiCounter,2} = ...
                        [nuclei_centroid_xy_i(1,1).Centroid(:,1), nuclei_centroid_xy_i(1,1).Centroid(:,2)];
                    
                    phi_new_masks_set{i,1}{nucleiCounter,1} = logical(clumpMask_j);
                    
                    GeoInitialGuessMaskSet_byImage{i,1}{nucleiCounter,1} = logical(clumpMask_j);
                    GeoInitialGuessMaskSet_byImage{i,1}{nucleiCounter,2} = ...
                        [nuclei_centroid_xy_i(1,1).Centroid(:,1), nuclei_centroid_xy_i(1,1).Centroid(:,2)];
                    
                    rough_DistMapSetByImage{i,1}{nucleiCounter,1} = -2 * ones(size(imSet{i,1}));
                    rough_DistMapSetByImage{i,1}{nucleiCounter,2} = ...
                        [nuclei_centroid_xy_i(1,1).Centroid(:,1), nuclei_centroid_xy_i(1,1).Centroid(:,2)];
                    
                    DistMapSetByImage_GeoCentroids{i,1}{nucleiCounter,1} = -2 * ones(size(imSet{i,1}));                    
                    DistMapSetByImage_GeoCentroids{i,1}{nucleiCounter,2} = ...
                        [nuclei_centroid_xy_i(1,1).Centroid(:,1), nuclei_centroid_xy_i(1,1).Centroid(:,2)];
                    
                    %+---------------------------------------+
                    %| Leave EMPTY element for Initial Guess |
                    %+---------------------------------------+
                    ClumpDistMapSetByImage{i,1}{j,1} = []; %-2 * ones(size(imSet{i,1}));                   
                    
                    nucleiCounter = nucleiCounter + 1;
                    continue;
                end

                %+------------------------------------------+
                %|           No nuclei cytoplasm            |
                %|  Should Find Local Minimal as Nuclei     |
                %| Possibly it's a false negative of nuclei |
                %+------------------------------------------+
                if insideNucleiNum == 0
                    continue;
                end

                %+------------------------------------------+
                %|      For Clump with Multiple Nucleus     |
                %|           Contains Overlapping           |
                %|      Need Initial Guess Computation!     |
                %+------------------------------------------+
                fprintf('\tCompute roughly distance map for each nuclei of Image %d, Clump %d\n', i, j);
                %+----------------------------------+
                %| Rough Initial Guess of Each Cell |
                %+----------------------------------+
                rough_DistMapSetByNuclei = compute_rad_distmap(clumpMask_j, nucleusMask_inside_Clump_j, beta); 

                rough_InitialGuessMaskSetByNuclei = cell(size(rough_DistMapSetByNuclei,1),2);  
                geo_InitialGuessMaskSetByNuclei = cell(size(rough_DistMapSetByNuclei,1),2);

                %+--------------------------+
                %| Covert rem to binary rem |
                %+--------------------------+
                for n  = 1:size(rough_DistMapSetByNuclei,1)
                    rough_InitialGuessMaskSetByNuclei{n,1} = false(size(rough_DistMapSetByNuclei{n,1}));
                    rough_InitialGuessMaskSetByNuclei{n,1}(rough_DistMapSetByNuclei{n,1} .* clumpMask_j ~= 0) = 1;
                    rough_InitialGuessMaskSetByNuclei{n,2} = rough_DistMapSetByNuclei{n,2};
                end

                %+----------------------------------+
                %|     Distance Map for A Clump     |
                %+----------------------------------+
                %| Include Process Distance Values  |
                %|      in Overlapping Regions      |
                %+----------------------------------+
                %| New Initial Guess for Each Cell  |
                %|      By Geometric Centroids      |
                %+----------------------------------+
                fprintf('\tCompute global distance map & distance map for each nuclei of image %d, clump %d by geometric centroids...\n', ...
                        i, j);
                [ClumpDistMapSetByImage{i,1}{j,1}, DistMapSetByNuclei_GeoCentroids] = ...
                                            BuildGlobalDistTerm(rough_InitialGuessMaskSetByNuclei, clumpMask_j, beta);
                

                %+--------------------------+
                %| Covert rem to binary rem |
                %+--------------------------+
                for n  = 1:size(DistMapSetByNuclei_GeoCentroids,1)
                    geo_InitialGuessMaskSetByNuclei{n,1} = false(size(DistMapSetByNuclei_GeoCentroids{n,1}));
                    geo_InitialGuessMaskSetByNuclei{n,1}(DistMapSetByNuclei_GeoCentroids{n,1} .* clumpMask_j ~= 0) = 1;
                    geo_InitialGuessMaskSetByNuclei{n,2} = DistMapSetByNuclei_GeoCentroids{n,2};    % Nuclei Centroid (x,y)
                    
                    %+----------------------------------+
                    %|      Remove tiny fragments       |
                    %+----------------------------------+
                    bwlabel_geo_InitialGuessMask = bwlabel(geo_InitialGuessMaskSetByNuclei{n,1},8);
                    if max(bwlabel_geo_InitialGuessMask(:)) > 1
                        for region_id = 1:max(bwlabel_geo_InitialGuessMask(:))
                            idx_label_regID = find(bwlabel_geo_InitialGuessMask == region_id);
                            if size(idx_label_regID,1) < min_TinyFragments_DistMap
                                geo_InitialGuessMaskSetByNuclei{n,1}(idx_label_regID) = 0;
                            end
                        end
                    end
                end

                t_distMap(i) = toc;                            
                fprintf('\tTime: %f\n', t_distMap(i));
                
                %+-------------------------+
                %|  Save Initial Results   |
                %+-------------------------+
                fprintf('\t\tSaving distance maps for image %d, clump %d...', i, j);

                for n = 1:size(rough_InitialGuessMaskSetByNuclei,1)                            
                    %+--------------------------------+
                    %|      Save Initial Guess        |
                    %+--------------------------------+
                    %| GeoInitialGuessMaskSet is used |
                    %|       in Next Section.         |
                    %+--------------------------------+
                    roughInitialGuessMaskSet_byImage{i,1}{nucleiCounter,1} = rough_InitialGuessMaskSetByNuclei{n,1};
                    roughInitialGuessMaskSet_byImage{i,1}{nucleiCounter,2} = rough_InitialGuessMaskSetByNuclei{n,2};
                    GeoInitialGuessMaskSet_byImage{i,1}{nucleiCounter,1} = geo_InitialGuessMaskSetByNuclei{n,1};
                    GeoInitialGuessMaskSet_byImage{i,1}{nucleiCounter,2} = geo_InitialGuessMaskSetByNuclei{n,2};
                    
                    %+--------------------------------+
                    %|   Save Initial Distance Map    |
                    %+--------------------------------+
                    %|          for Each Cell         |
                    %|               or               |
                    %|        Global of a Clump       |
                    %+--------------------------------+
                    %| DistMapSetByImage_GeoCentroids |
                    %|    is used in Next Section.    |
                    %+--------------------------------+
                    rough_DistMapSetByImage{i,1}{nucleiCounter,1} = rough_DistMapSetByNuclei{n,1};
                    rough_DistMapSetByImage{i,1}{nucleiCounter,2} = rough_DistMapSetByNuclei{n,2};
                    DistMapSetByImage_GeoCentroids{i,1}{nucleiCounter,1} = DistMapSetByNuclei_GeoCentroids{n,1};
                    DistMapSetByImage_GeoCentroids{i,1}{nucleiCounter,2} = DistMapSetByNuclei_GeoCentroids{n,2};

                    nucleiCounter = nucleiCounter + 1;
                end

                fprintf('Saved!\n');                        
            end
            
            pause(0.1);
        end

        %+-------------------------------+
        %| Save section variable to file |
        %+-------------------------------+
        save(strcat(OUTPUT_PATH, storageInitial, '/initialPhis_beta_', num2str(beta), '.mat'), ...
            'ClumpDistMapSetByImage', 'GeoInitialGuessMaskSet_byImage', 'DistMapSetByImage_GeoCentroids', 't_distMap', ...
            'roughInitialGuessMaskSet_byImage', 'rough_DistMapSetByImage');        
        pause(60);

        fprintf('Clearing results for beta = %f...\n', beta);
        clear ClumpDistMapSetByImage GeoInitialGuessMaskSet_byImage DistMapSetByImage_GeoCentroids t_distMap ...
              roughInitialGuessMaskSet_byImage rough_DistMapSetByImage;
        
        pause(60);
    end
end


%% 
%===============================================
% Overlapping segmentation by distance map term
%===============================================
try
    load('xxxx');
    fprintf('Step - 6 LSF evolution by initial guess.........done!\n');
catch
    fprintf('Step - 6 LSF evolution by initial guess...\n');
    
    for para_id_BETA = 1:size(beta_logistic_set,1)
        beta = beta_logistic_set(para_id_BETA);
        
        load(strcat(OUTPUT_PATH, storageInitial, 'initialPhis_beta_', num2str(beta), '.mat'), ...
            'ClumpDistMapSetByImage', 'GeoInitialGuessMaskSet_byImage')
           
        fprintf('+--------------------------+\n');
        fprintf('|     Beta = %f      |\n', beta);
        fprintf('+--------------------------+\n');
        
         for kappa_setID = 1:length(kappa_set)
             kappa_LSF = kappa_set(kappa_setID);

            for chaID = 1:length(chi_set)
                chi_LSF = chi_set(chaID);

                fprintf('Kappa = %f, Chi = %f\n', kappa_LSF, chi_LSF);
                    %==================================
                    %         Global variables
                    %==================================
                    InitialGuessMaskSet_inComputing = GeoInitialGuessMaskSet_byImage;
                    ClumpDistMapsSet_inComputing = ClumpDistMapSetByImage;

                    %==================================
                    %  Mat for different loop numbers
                    %==================================               
                    LSF_5 = cell(imNum,1);
                    LSF_4 = cell(imNum,1);
                    LSF_3 = cell(imNum,1);
                    LSF_2 = cell(imNum,1);
                    LSF_1 = cell(imNum,1);

                    ClumpDistMapsSet_inComputing_1 = cell(imNum,1);
                    ClumpDistMapsSet_inComputing_2 = cell(imNum,1);
                    ClumpDistMapsSet_inComputing_3 = cell(imNum,1);
                    ClumpDistMapsSet_inComputing_4 = cell(imNum,1);

                    for i = 1:imNum
                        tic;
                        fprintf('\tLSF Evolution, image %d\n', i);
                        num_loop = 1;
                        while(num_loop <= loop)
                            im_i = imCytoSet{i,1};
                            %+------------------------------------------+
                            %| Initial Guesses for ALL Cells in Image i |
                            %+------------------------------------------+
                            InitialGuessMaskSet_i = InitialGuessMaskSet_inComputing{i,1};

                            fprintf('\tLoop info: Image %d, %d-th loop.\n', i, num_loop);

                            %+------------------------------+
                            %| ALL Clumps Masks for Image i |
                            %+------------------------------+
                            ClumpsSetMask_i = SceneCytoClumpMaskSet{i,1};   % clumps mask
                            clumpsSetStats_i = regionprops(ClumpsSetMask_i, 'PixelIdxList');  % pixels lists of each clump

                            %==================================
                            % For each cell and its neighbours
                            %==================================
                            for j = 1:size(InitialGuessMaskSet_i,1)
                                phi_j = InitialGuessMaskSet_i{j,1};
                                phi_j_Idx = find(phi_j == 1);

                                %+----------------------------------+
                                %| If phi_1 MISSED in LSF Evolution | 
                                %+----------------------------------+
                                if isempty(phi_j_Idx)
                                    continue;
                                end

                                %========================================
                                % Find the clump contains this cytoplasm
                                %========================================
                                clumpContainPhi_j = false(size(im_i));
                                I = im_i;
                                for k = 1:size(clumpsSetStats_i,1)
                                    clumpPixelsIdx_k = clumpsSetStats_i(k,1).PixelIdxList;

                                    if ~isempty(intersect(phi_j_Idx, clumpPixelsIdx_k))
                                        clumpContainPhi_j(clumpPixelsIdx_k) = 1;
                                        I( clumpContainPhi_j ~= 1 ) = 0;
                                        break;
                                    end
                                end

                                %+---------------------------------+
                                %|     Find phi_1's neighbours     |
                                %+---------------------------------+                        
                                neighborsOfPhi_j = cell(1,1);                        
                                neighborNum = 1;
                                for k = 1:size(InitialGuessMaskSet_i,1)
                                    if k == j                    
                                        continue;
                                    end
                                    phi_k_TempImg = InitialGuessMaskSet_i{k,1};
                                    phi_k_PixelsIdx = find(phi_k_TempImg == 1);

                                    %+----------------------------------+
                                    %| If phi_2 MISSED in LSF Evolution | 
                                    %+----------------------------------+
                                    if isempty(phi_k_PixelsIdx)
                                        continue;
                                    end

                                    if ~isempty(intersect(phi_j_Idx, phi_k_PixelsIdx))
                                        neighborsOfPhi_j{neighborNum,1} = phi_k_TempImg;
                                        neighborNum = neighborNum + 1;
                                    else
                                        % check if boundaries are close enough
                                        bwboundaries_j = bwboundaries(phi_j);
                                        bwboundaries_k = bwboundaries(phi_k_TempImg);

                                        boundary_j = bwboundaries_j{1,1};
                                        boundary_k = bwboundaries_k{1,1};
                                        distsListPts = zeros(size(boundary_j,1) * size(boundary_k,1),1);
                                        numDist = 1;
                                        for p = 1:size(boundary_j,1)
                                            for q = 1:size(boundary_k,1)
                                                distsListPts(numDist,1) = sqrt((boundary_j(p,1) - boundary_k(q,1))^2 ...
                                                                                + (boundary_j(p,2) - boundary_k(q,2))^2);
                                                numDist = numDist + 1;
                                            end
                                        end

                                        minDistsPt = min(distsListPts(:));
                                        if minDistsPt <= min_PtDistance
                                            neighborsOfPhi_j{neighborNum,1} = phi_k_TempImg;
                                            neighborNum = neighborNum + 1;
                                        end
                                    end
                                end

                                %+--------------------------------+
                                %|     Single Cell in a Clump     |
                                %| Push it into the Cell Mask Set |
                                %+--------------------------------+
                                if length(neighborsOfPhi_j) == 1 && isempty(neighborsOfPhi_j{1,1}) == 1
                                    InitialGuessMaskSet_inComputing{i,1}{j,1} = phi_j;
                                    continue;
                                end

                                %==================================
                                % Find the pair dist map of phi_1
                                %==================================
                                for k = 1:size(ClumpDistMapsSet_inComputing{i,1},1);
                                    ClumpDistMap_k = ClumpDistMapsSet_inComputing{i,1}{k,1};
                                    if isempty(ClumpDistMap_k)
                                        continue;
                                    end
                                    if max(ClumpDistMapsSet_inComputing{i,1}{k,1}(:)) == 0
                                        continue;
                                    end
                                    idx_map = find(ClumpDistMap_k ~= 1);
                                    if ~isempty(intersect(idx_map, phi_j_Idx))
                                        ClumpDistMap_ContainPhi_j = ClumpDistMap_k;
                                        break;
                                    end
                                end

                                %+----------------------------+ 
                                %|      *LSF EVOLUTION*       |
                                %|(cell pair: phi_j and phi_k)|
                                %+----------------------------+
                                phi_LSF_j = phi_j;
                                for k = 1:length(neighborsOfPhi_j)
                                    phi_LSF_k = neighborsOfPhi_j{k,1};
                                    if islogical(phi_LSF_j) == 0
                                        phi_LSF_j = ~im2bw(phi_LSF_j);
                                    end

                                    %+------------------------------+
                                    %|       LSF: Update phi_j      |
                                    %+------------------------------+
                                    fprintf('\tCell: %d, Neighbour: %d\n',  j, k)
                                    [phi_LSF_j] = ...
                                                distLevelSet(I, phi_LSF_j, phi_LSF_k, ClumpDistMap_ContainPhi_j, ...
                                                iter_in_extent, iter_out_extent, ...
                                                kappa_LSF, chi_LSF);

                                    %+------------------------------+
                                    %| phi_LSF_j: double to logical |
                                    %+------------------------------+
                                    phi_j_logical_tmp = true(size(phi_LSF_j));
                                    phi_j_logical_tmp(phi_LSF_j >= 0) = 0;
                                    phi_LSF_j = phi_j_logical_tmp;

                                    %+------------------------------+
                                    %|  Filling holes in phi_LSF_j  |
                                    %+------------------------------+
                                    InitialGuessMaskSet_inComputing{i,1}{j,1} = imfill(phi_LSF_j, 'holes');

                                    %+---------------------------------------------+
                                    %| Remove smaller regions in phi_j (fragments) |
                                    %+---------------------------------------------+
                                    phi_j_labels = bwlabel(InitialGuessMaskSet_inComputing{i,1}{j,1}, 8);
                                    if max(phi_j_labels(:)) > 1
                                        phi_j_labels_tmp = false(size(InitialGuessMaskSet_inComputing{i,1}{j,1}));
                                        regionSize_phi_j = zeros(1,1);
                                        for label_id = 1:max(phi_j_labels(:))
                                            regionSize_phi_j(label_id,1) = size(find(phi_j_labels == label_id),1);
                                        end
                                        maxRegion_phi_j_id = find(regionSize_phi_j == max(regionSize_phi_j(:)));
                                        idx_max_label = find(phi_j_labels == maxRegion_phi_j_id);
                                        phi_j_labels_tmp(idx_max_label) = 1;
                                        InitialGuessMaskSet_inComputing{i,1}{j,1} = phi_j_labels_tmp;
                                    end
                                end
                            end

                    %+---------------------------------------------------+
                    %         [ABOVE] LSF Evolution for each cell
                    %+---------------------------------------------------+
                    %+---------------------------------------------------+
                    %          [BELOW] Update Clump Distance Map
                    %+---------------------------------------------------+

                            idx_visitedCell = zeros(1,1);   % cells in clumps and visited
                            ClumpsSetMask_i = SceneCytoClumpMaskSet{i,1};   % clumps mask
                            clumpsSetStats_i = regionprops(ClumpsSetMask_i, 'PixelIdxList');  % pixels lists of each clump

                            %================================
                            % For each cell and its partners
                            %================================
                            for j = 1:size(InitialGuessMaskSet_inComputing{i,1},1)

                                %================================
                                %    Cell j visited, continue
                                %================================
                                if ~isempty(intersect(j, idx_visitedCell))
                                    continue;
                                end

                                phi_j = InitialGuessMaskSet_inComputing{i,1}{j,1};
                                phi_j_Idx = find(phi_j == 1);

                                %========================================
                                % Find the clump contains this cytoplasm
                                %========================================
                                clumpContainPhi_j = zeros(size(im_i));
                                for k = 1:length(clumpsSetStats_i)
                                    clumpPixelsIdx_k = clumpsSetStats_i(k,1).PixelIdxList;
                                    isNotInteract = isempty(intersect(phi_j_Idx, clumpPixelsIdx_k));

                                    if isNotInteract == 0
                                        clumpID = k;
                                        clumpContainPhi_j(clumpPixelsIdx_k) = 1;
                                        clumpContainPhi_j = logical(clumpContainPhi_j);
                                        break;
                                    end
                                end

                                %===============================================
                                % Find phi's partners inside 'clumpContainPhi1'
                                %===============================================
                                phiClumpPartners = cell(1,1);
                                numPartner = 1;
                                for k = 1:size(InitialGuessMaskSet_inComputing{i,1},1)
                                    if k == j
                                        phiClumpPartners{numPartner,1} = InitialGuessMaskSet_inComputing{i,1}{k,1};
                                        phiClumpPartners{numPartner,2} = InitialGuessMaskSet_inComputing{i,1}{k,2};
                                        numPartner = numPartner + 1;
                                        idx_visitedCell(end + 1,1) = k;
                                        continue;
                                    end

                                    clump_idx = find(clumpContainPhi_j == 1);
                                    phi_partner_idx = find(InitialGuessMaskSet_inComputing{i,1}{k,1} == 1);
                                    if ~isempty(intersect(clump_idx,phi_partner_idx))
                                        phiClumpPartners{numPartner,1} = InitialGuessMaskSet_inComputing{i,1}{k,1};
                                        phiClumpPartners{numPartner,2} = InitialGuessMaskSet_inComputing{i,1}{k,2};
                                        numPartner = numPartner + 1;
                                        idx_visitedCell(end + 1,1) = k;
                                    end
                                end

                                %+-------------------------------+
                                %|       For single-cell,        |
                                %| push it into the clumps stack |
                                %+-------------------------------+
                                if size(phiClumpPartners,1) <= 1
                                    continue;
                                end                            

                                %+---------------------------------+
                                %|      Update pair dist map       |
                                %+---------------------------------+
                                fprintf('\tUpdating global distance map for Clump %d through Cell %d...\n', clumpID, j);
                                ClumpDistMap_ContainPhi_j = BuildGlobalDistTerm(phiClumpPartners, clumpContainPhi_j, beta);

                                ClumpDistMapsSet_inComputing{i,1}{clumpID,1} = ClumpDistMap_ContainPhi_j;
                            end

                            LSF_5{i,1} = InitialGuessMaskSet_inComputing{i,1};                                

                            t_LSFSeg(i) = toc;

                            %+--------------------------------------+
                            %| Assign LSF result to stage variables |
                            %+--------------------------------------+
                            switch num_loop
                                case 1
                                    LSF_1{i,1} = LSF_5{i,1};
                                    ClumpDistMapsSet_inComputing_1{i,1} = ClumpDistMapsSet_inComputing{i,1};
                                    t_LSFSeg_1(i) = t_LSFSeg(i);
                                case 2
                                    LSF_2{i,1} = LSF_5{i,1};
                                    ClumpDistMapsSet_inComputing_2{i,1} = ClumpDistMapsSet_inComputing{i,1};
                                    t_LSFSeg_2(i) = t_LSFSeg(i);
                                case 3
                                    LSF_3{i,1} = LSF_5{i,1};
                                    ClumpDistMapsSet_inComputing_3{i,1} = ClumpDistMapsSet_inComputing{i,1};
                                    t_LSFSeg_3(i) = t_LSFSeg(i);
                                case 4
                                    LSF_4{i,1} = LSF_5{i,1};
                                    ClumpDistMapsSet_inComputing_4{i,1} = ClumpDistMapsSet_inComputing{i,1};
                                    t_LSFSeg_4(i) = t_LSFSeg(i);                           
                            end

                            %=====================================
                            %           Update num_loops
                            %=====================================
                            num_loop = num_loop + 1;
                        end

                        fprintf('\tLSF Time: %f\n', t_LSFSeg(i));
                    end
                    
                    %==================================
                    %  Save section variable to file
                    %==================================
                    save(strcat(OUTPUT_PATH, storageDist, 'LSF5/LSF_5_beta_', num2str(beta), ...
                        '_kappa_', num2str(kappa_LSF), '_chi_', num2str(chi_LSF), '_iterIn_', num2str(iter_in_extent), '_iterOut_', num2str(iter_out_extent), '.mat'),...
                                'LSF_5', 'ClumpDistMapsSet_inComputing', 't_LSFSeg');

                    fprintf('\t\tSaving LSF5 ...\n');
                    pause(60);

                    save(strcat(OUTPUT_PATH, storageDist, 'LSF4/LSF_4_beta_', num2str(beta), ...
                        '_kappa_', num2str(kappa_LSF), '_chi_', num2str(chi_LSF), '_iterIn_', num2str(iter_in_extent), '_iterOut_', num2str(iter_out_extent), '.mat'),...
                                'LSF_4', 'ClumpDistMapsSet_inComputing_4', 't_LSFSeg_4');

                    fprintf('\t\tSaving LSF4 ...\n');
                    pause(60);

                    save(strcat(OUTPUT_PATH, storageDist, 'LSF3/LSF_3_beta_', num2str(beta), ...
                        '_kappa_', num2str(kappa_LSF), '_chi_', num2str(chi_LSF), '_iterIn_', num2str(iter_in_extent), '_iterOut_', num2str(iter_out_extent), '.mat'),...
                                'LSF_3', 'ClumpDistMapsSet_inComputing_3', 't_LSFSeg_3');

                    fprintf('\t\tSaving LSF3 ...\n');
                    pause(60);

                    save(strcat(OUTPUT_PATH, storageDist, 'LSF2/LSF_2_beta_', num2str(beta), ...
                        '_kappa_', num2str(kappa_LSF), '_chi_', num2str(chi_LSF), '_iterIn_', num2str(iter_in_extent), '_iterOut_', num2str(iter_out_extent), '.mat'),...
                                'LSF_2', 'ClumpDistMapsSet_inComputing_2', 't_LSFSeg_2');

                    fprintf('\t\tSaving LSF2 ...\n');
                    pause(60);

                    save(strcat(OUTPUT_PATH, storageDist, 'LSF1/LSF_1_beta_', num2str(beta), ...
                        '_kappa_', num2str(kappa_LSF), '_chi_', num2str(chi_LSF), '_iterIn_', num2str(iter_in_extent), '_iterOut_', num2str(iter_out_extent), '.mat'),...
                                'LSF_1', 'ClumpDistMapsSet_inComputing_1', 't_LSFSeg_1');

                    fprintf('\t\tSaving LSF1 ...\n');
                    pause(60);
            end
        end
        
        fprintf('Clearing results for beta = %f, kappa = %f, chi = %f ...\n', beta, kappa_LSF, chi_LSF);
				
        clear ClumpDistMapSetByImage GeoInitialGuessMaskSet_byImage
        pause(60);
    end
end
end
