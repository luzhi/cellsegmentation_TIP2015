function [ CentroidBoundaryPtsPairSet, MinimalPartitionValue, XYPtsOnLinkLines ] = ...
            ComputeCentroidBoundaryPairs( ClumpMask, nucleusMaskInsideClump )
% Compute "Centroid-Boundary Points Pairs" for ALL cells in the clump.
%
%   "Boundary" here means either the clump boundary or the boundary of a
%   group of Internal Cells.
%
%   *Note: This is a "recursive" function.


    %+-------------------+
    %| Nucleus Centroids |
    %+-------------------+
    nucleusRegions = regionprops(nucleusMaskInsideClump, 'Centroid', 'PixelIdxList');
    CentroidsXY = zeros(length(nucleusRegions),2);
    for i = 1:size(nucleusRegions,1)
        CentroidsXY(i,1) = nucleusRegions(i,1).Centroid(:,1);      % column - x
        CentroidsXY(i,2) = nucleusRegions(i,1).Centroid(:,2);      % row - y    
    end

    %+-----------------------+
    %| Clump Boundary Pixels |
    %+-----------------------+
    objBoundary = bwboundaries(ClumpMask);
    BoundaryPixels(:,1) = objBoundary{1,1}(:,2);     % x-axis <--> column
    BoundaryPixels(:,2) = objBoundary{1,1}(:,1);     % y-axis <--> row

    %+----------------------------------------------+
    %| Linking Lines without those in Lapped Region |
    %|        (Centroids, Boundary_Points)          |
    %+----------------------------------------------+
    %|       By "Nearest Centroid" Principle        |
    %+----------------------------------------------+
    [CentroidBoundaryPtsPairSet, InternalCentroidIndicator, MinimalPartitionValue, XYPtsOnLinkLines] ...
                                                = SelectNearestCentroid(CentroidsXY, BoundaryPixels, ClumpMask);
    
    %+------------------------------------+
    %|      Delete Duplicate Points       |
    %+------------------------------------+
    RemoveIDList = [];
    numRemoveID = 1;
    for i = 1:(size(CentroidBoundaryPtsPairSet,1) - 1)
        if ~isempty(RemoveIDList)
            if ismember(i, RemoveIDList)
                continue;
            end
        end
        
        Pt_i(1,1:2) = CentroidBoundaryPtsPairSet{i,1};
        Pt_i(1,3:4) = CentroidBoundaryPtsPairSet{i,2};
        for j = (i + 1):size(CentroidBoundaryPtsPairSet,1)        
%             if i == j
%                 continue;
%             end

            if ~isempty(RemoveIDList)
                if ismember(j, RemoveIDList)
                    continue;
                end
            end

            Pt_j(1,1:2) = CentroidBoundaryPtsPairSet{j,1};
            Pt_j(1,3:4) = CentroidBoundaryPtsPairSet{j,2};
            if isequal(Pt_i, Pt_j)
               RemoveIDList(numRemoveID,1) = j;
               numRemoveID = numRemoveID + 1;
               isSingle = 0;
            end
        end
    end

    CentroidBoundaryPtsPairSet(RemoveIDList(:),:) = [];
    XYPtsOnLinkLines(RemoveIDList(:),:) = [];
    
    %+-----------------------------------+
    %|  Check if only ONE linking line?  |
    %+-----------------------------------+
    [tmp_CBPtsPariSet, ~] =  GroupCBPtsPairsByCentroids(CentroidBoundaryPtsPairSet, XYPtsOnLinkLines, CentroidsXY);
    for i = 1:size(tmp_CBPtsPariSet,1)
        if size(tmp_CBPtsPariSet{i,1},1) <= 1
            InternalCentroidIndicator(i,1) = 1;
            centroidxy_i = CentroidsXY(i,1:2);
            idx_tobe_removed = [];
            num_tobe_removed = 1;
            for j = 1:size(CentroidBoundaryPtsPairSet,1)
                centroidxy_j = CentroidBoundaryPtsPairSet{j,1};
                if isequal(centroidxy_i, centroidxy_j)
                    idx_tobe_removed(num_tobe_removed,1) = j;
                end
            end
            
            CentroidBoundaryPtsPairSet(idx_tobe_removed(:),:) = [];
        end
    end
    
    %+----------------------------------------------+
    %|        Linking Lines in Lapped Region        |
    %+----------------------------------------------+
    %|             By "Vector Rotation"             |
    %+----------------------------------------------+
    [CentroidBoundaryPtsPairSet, XYPtsOnLinkLines] = ...
            RotateVectorCrossLappedRegion(CentroidBoundaryPtsPairSet, XYPtsOnLinkLines);
        
    %+-------------------------------------+
    %|     Case 1: Has Internal Cells      |
    %+-------------------------------------+
    if ismember(1, InternalCentroidIndicator)
        %+------------------------------+
        %|    Internal Nucleus Mask     |
        %+------------------------------+
        idx_InternalCells = find(InternalCentroidIndicator == 1);
        nucleusMask_Internal = false(size(ClumpMask));
        for i = 1:size(idx_InternalCells,1)
            idx_InternalCentroid_PixelList = ...
                    nucleusRegions(idx_InternalCells(i,1),1).PixelIdxList(:);
            nucleusMask_Internal(idx_InternalCentroid_PixelList) = 1;
        end
        
        %+---------------------------------+
        %| Distance among Internal Nucleus |
        %+---------------------------------+
        %|     (Get the minimal one!)      |
        %+---------------------------------+
        dist_among_internal_nucleus = zeros(1,3);
        numDist = 1;
        %+--------------------------------------+
        %|  If there is only ONE internal cell  |
        %+--------------------------------------+
        %| Distance = internal to nearest outer |
        %+--------------------------------------+        
        if size(idx_InternalCells,1) == 1
            idx_OuterCells = find(InternalCentroidIndicator == 0);
            centroidsxy_inner(1,1) = CentroidsXY(idx_InternalCells(1,1),1);
            centroidsxy_inner(1,2) = CentroidsXY(idx_InternalCells(1,1),2);
            dist_among_internal_outer = zeros(1,1);
            for i = 1:size(idx_OuterCells,1)
                centroidsxy_outer_i(1,1) = CentroidsXY(idx_OuterCells(i,1),1);
                centroidsxy_outer_i(1,2) = CentroidsXY(idx_OuterCells(i,1),2);
                
                dist_among_internal_outer(i,1) = pdist([centroidsxy_inner; centroidsxy_outer_i], 'euclidean');
            end
            
            id_nearest_outer_cell = find(dist_among_internal_outer == min(dist_among_internal_outer(:)));
            dist_among_internal_nucleus(numDist,1) = min(dist_among_internal_outer(:));
            dist_among_internal_nucleus(numDist,2) = idx_InternalCells(1,1);
            dist_among_internal_nucleus(numDist,3) = idx_OuterCells(id_nearest_outer_cell(1,1),1);
        else
            %+----------------------------------------+
            %|   If there is MULTIPLE internal cell   |
            %+----------------------------------------+
            %| Distance = internal to nearest another |
            %+----------------------------------------+
            for i = 1:(size(idx_InternalCells,1) - 1)
                centroidsxy_i(1,1) = CentroidsXY(idx_InternalCells(i,1),1);
                centroidsxy_i(1,2) = CentroidsXY(idx_InternalCells(i,1),2);

                for j = (i + 1):size(idx_InternalCells,1)
                    centroidsxy_j(1,1) = CentroidsXY(idx_InternalCells(j,1),1);
                    centroidsxy_j(1,2) = CentroidsXY(idx_InternalCells(j,1),2);

                    dist_among_internal_nucleus(numDist,1) = pdist([centroidsxy_i; centroidsxy_j], 'euclidean');
                    dist_among_internal_nucleus(numDist,2) = idx_InternalCells(i,1);
                    dist_among_internal_nucleus(numDist,3) = idx_InternalCells(j,1);
                    numDist = numDist + 1;
                end
            end        
            dist_among_internal_nucleus = [dist_among_internal_nucleus; ...
                                           dist_among_internal_nucleus(:,1), dist_among_internal_nucleus(:,3), dist_among_internal_nucleus(:,2)];
        end
        
        %+---------------------------------------+
        %| Get "boundary points" for each Nuclei |
        %+---------------------------------------+
        %|  Generate points on the virtual line. |
        %+---------------------------------------+
        ESTIMATED_RATIO = 0.70;
        for i = 1:size(idx_InternalCells,1)
            %+---------------------------------------+
            %| Get "boundary points" for each Nuclei |
            %+---------------------------------------+
            id_internal_cell_i = idx_InternalCells(i,1);
            idx_i_dist = find(dist_among_internal_nucleus(:,2) == id_internal_cell_i);
            dist_internal_cell_i = dist_among_internal_nucleus(idx_i_dist(:),:);
            
            if size(dist_internal_cell_i,1) == 1
                NumPointsOnLine = floor(ESTIMATED_RATIO * dist_internal_cell_i(1,1));
                centroidsxy_j(1,1) = CentroidsXY(dist_internal_cell_i(1,3),1);
                centroidsxy_j(1,2) = CentroidsXY(dist_internal_cell_i(1,3),2);
            end
            
            if size(dist_internal_cell_i,1) > 1
                NumPointsOnLine = floor(ESTIMATED_RATIO * min(dist_internal_cell_i(:,1)));
                idx_min_j = find(dist_internal_cell_i == min(dist_internal_cell_i(:,1)));
                centroidsxy_j(1,1) = CentroidsXY(dist_internal_cell_i(idx_min_j(1,1),3),1);
                centroidsxy_j(1,2) = CentroidsXY(dist_internal_cell_i(idx_min_j(1,1),3),2);
            end            
            
            MinimalPartitionValue = min(MinimalPartitionValue, NumPointsOnLine);            
            
            %+---------------------------------------+
            %|  Generate points on the virtual line. |
            %+---------------------------------------+
            centroidsxy_i(1,1) = CentroidsXY(idx_InternalCells(i,1),1);
            centroidsxy_i(1,2) = CentroidsXY(idx_InternalCells(i,1),2);
            
            Length_i = pdist([centroidsxy_i; centroidsxy_j], 'euclidean');
            
            CBPtsPairSetofaCell = cell(1,2);
            numCBPtsOfaCell = 1;
            
            Vec_i = [centroidsxy_j(1,1) - centroidsxy_i(1,1); ...
                     centroidsxy_j(1,2) - centroidsxy_i(1,2)];
            iVec_i = Vec_i / sqrt(sum(Vec_i .^2, 1));
            
            num_rot_step = 359;
            
            for j = 1:num_rot_step
                M_rot_step = [cos(j * (pi / 180) *1), -sin(j * (pi / 180) *1); ...
                                     sin(j * (pi / 180) *1),  cos(j * (pi / 180) *1)];
                iVec_rot_j = M_rot_step * iVec_i;
                
                Vec_rot_j =  (ESTIMATED_RATIO * Length_i) .* iVec_rot_j;
                Vec_rot_j = Vec_rot_j';
                
                IS_NaN = isnan(Vec_rot_j);
                IS_Inf = isinf(Vec_rot_j);
                if ~ismember(1,IS_NaN) && ~ismember(1,IS_Inf)
                    CBPtsPairSetofaCell{numCBPtsOfaCell,1} = centroidsxy_i(1,1:2);
                    CBPtsPairSetofaCell{numCBPtsOfaCell,2} = ...
                        [centroidsxy_i(1,1) + Vec_rot_j(1,1), centroidsxy_i(1,2) + Vec_rot_j(1,2)];                
                    numCBPtsOfaCell = numCBPtsOfaCell + 1;
                end
            end
            CentroidBoundaryPtsPairSet = [CentroidBoundaryPtsPairSet; CBPtsPairSetofaCell];
            
            %+--------------------------------+
            %|          Get Results           |
            %+--------------------------------+
            for j = 1:size(CBPtsPairSetofaCell,1)
%                 j
%                                          CBPtsPairSetofaCell{j,1}
%                                          CBPtsPairSetofaCell{j,2}
                                         
                XYPtsOnLinkLines = [XYPtsOnLinkLines; ...
                    GetInterpolationPointsXY(CBPtsPairSetofaCell{j,1}, CBPtsPairSetofaCell{j,2}, ...
                                             NumPointsOnLine)];                                        
            end
        end
    end
    
    if MinimalPartitionValue == 0
        MinimalPartitionValue = 1;
    end
end

