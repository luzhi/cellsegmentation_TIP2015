function varargout = BuildGlobalDistTerm(InitialGuessMaskSet_byNuclei, ClumpMask, beta)
% Build distance map for a pair of overlapping cells.

    numCells = size(InitialGuessMaskSet_byNuclei,1);
    distTerm = zeros(size(InitialGuessMaskSet_byNuclei{1,1}));
    phis_overlapping = [];
    numOverRegion = 1;

    phis_initial = InitialGuessMaskSet_byNuclei;

    %====================================
    % Compute new centroids for each phi
    %====================================
    cens = zeros(numCells,2);
    for i = 1:size(phis_initial,1)
        phi_i = phis_initial{i,1};
        phiStats = regionprops(phi_i, 'Centroid', 'PixelList');
        if size(phiStats,1) > 1
            idx_length_labels = zeros(size(phiStats,1),1);
            for j = 1:size(phiStats,1)
                idx_length_labels(j,1) = size(phiStats(j,1).PixelList,1);
            end
            idx_max_length_labels = find(idx_length_labels == max(idx_length_labels));
            cens(i,1) = round(phiStats(idx_max_length_labels(1,1),1).Centroid(:,1));   % x - column
            cens(i,2) = round(phiStats(idx_max_length_labels(1,1),1).Centroid(:,2));   % y - row
            cens(i,3:4) = phis_initial{i,2};  % (x',y') of nuclei centroid (biology property)
        else            
            cens(i,1) = round(phiStats(1,1).Centroid(:,1));    % x - column
            cens(i,2) = round(phiStats(1,1).Centroid(:,2));    % y - row
            cens(i,3:4) = phis_initial{i,2};  % (x',y') of nuclei centroid (biology property)
        end
    end

    %+----------------------------------+
    %|    Compute new distance maps     |
    %| for each phi by its new centroid |
    %+----------------------------------+
    CentroidsMask = false(size(ClumpMask));
    for i = 1:size(cens,1)
%         CentroidsMask(round(cens(i,2)), round(cens(i,1))) = 1;
        CentroidsMask(cens(i,2), cens(i,1)) = 1;
    end
    phis_geo = compute_rad_distmap(ClumpMask, CentroidsMask, beta);
    
    %+----------------------------------+
    %|  Mapping Geo Centroid to Nuclei  |
    %+----------------------------------+
    %|        For Evaluation Code       |
    %+----------------------------------+
    assert(size(phis_geo,1) == size(cens,1), ...
        '[BuildGlobalDistTerm.m] phis_geo should have same size of cens');
    for i = 1:size(phis_geo,1)
        for j = 1:size(cens,1)
            if isequal(phis_geo{i,2}, cens(j,1:2))
                phis_geo{i,2} = cens(j,3:4);
                break;
            end
        end
    end
    
    %===========================
    %   Convert rem to binary
    %===========================
    phis_bw = cell(size(phis_geo,1),1);
    for i = 1:size(phis_geo,1)
        phis_bw{i,1} = false(size(phis_geo{i,1}));
        phis_bw{i,1}( phis_geo{i,1} .* ClumpMask ~= 0) = 1;
    end
    
    %===================================
    %     Find overlapping regions
    %===================================
    for i = 1:(size(phis_bw, 1) - 1)
        idx_i = find(phis_bw{i,1} == 1);
        for j = (i + 1): size(phis_bw, 1)
            idx_j = find(phis_bw{j,1} == 1);
            if ~isempty(intersect(idx_i, idx_j))
                idx_comm_ij = intersect(idx_i, idx_j);
                if size(idx_comm_ij,1) <= 5
                    continue;
                end
                phi_over = zeros(size(phis_bw{1,1}));
                phi_over(idx_comm_ij) = 1;
                phis_overlapping{numOverRegion,1} = logical(phi_over);
                numOverRegion = numOverRegion + 1;
            end
        end
    end

    %====================================================
    % Compute distance values in each overlapping reigon
    %====================================================
    if size(phis_overlapping,1) > 0
        for i = 1:size(phis_geo, 1)
            distTerm = distTerm + phis_geo{i,1};
        end
        
        %================================================
        % Choose values of pixels in overlapping reigons
        %================================================
        for i = 1:size(phis_overlapping,1)
            idx_comm_ith = phis_overlapping{i,1};
            maxValuesAtRegion = zeros(size(phis_overlapping{i,1}));
            for j = 1:size(phis_geo,1)
                maxValuesAtRegion(idx_comm_ith == 1) = ...
                    max(maxValuesAtRegion(idx_comm_ith == 1), ...
                        phis_geo{j,1}(idx_comm_ith == 1));
            end
            
            distTerm(idx_comm_ith == 1) = maxValuesAtRegion(idx_comm_ith == 1);
        end
    else
        for i = 1:size(phis_geo, 1)
            distTerm = distTerm + phis_geo{i,1};
        end
    end

    %==================================
    %   Remove values outside clump
    %==================================
    distTerm = distTerm .* ClumpMask;
    
    %====================================
    % Assign lagest values outside clump
    %====================================
    distTerm(ClumpMask ~= 1) = 1;
    
    varargout{1,1} = distTerm;
    varargout{2,1} = phis_geo;    
end
