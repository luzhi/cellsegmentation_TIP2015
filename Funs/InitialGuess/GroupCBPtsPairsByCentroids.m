function [CentroidBoundaryPtsPairSetByNuclei, XYPtsOnLinkLinesByNuclei] = ...
            GroupCBPtsPairsByCentroids(CentroidBoundaryPtsPairSet, XYPtsOnLinkLines, CentroidsXY)
% Rank the "Centroid-Boundary Points Pairs" by Centroids.
%
%   Centroid-Boundary Points Pairs with same centroid is grouped together.


    CentroidBoundaryPtsPairSetByNuclei = cell(size(CentroidsXY,1),1);
    XYPtsOnLinkLinesByNuclei = cell(size(CentroidsXY,1),1);
    
    %+-------------------------------------------------+
    %|               List of Centroids                 |
    %| in unranked Centroid-Boundary Points Pairs List |
    %+-------------------------------------------------+
    CBPtsPairsCentroids = zeros(size(CentroidBoundaryPtsPairSet,1),2);
    for i = 1:size(CentroidBoundaryPtsPairSet,1)
        CBPtsPairsCentroids(i,:) = CentroidBoundaryPtsPairSet{i,1};
    end
    
    %+-------------------+
    %| Rank by Centroids |
    %+-------------------+
    for i = 1:size(CentroidsXY,1)
        CentroidXY_i = CentroidsXY(i,:);
        membershipList = ismember(CBPtsPairsCentroids, CentroidXY_i, 'rows'); %ismember(CBPtsPairsCentroids, CentroidXY_i);
        idx_CentroidXY_i = find(membershipList(:,1) == 1);
        
        tmp_CentroidBoundaryPtsPairSetByNuclei = CentroidBoundaryPtsPairSet(idx_CentroidXY_i,:);
        CentroidBoundaryPtsPairSetByNuclei{i,1} = tmp_CentroidBoundaryPtsPairSetByNuclei;
        
        tmp_XYPtsOnLinkLines = XYPtsOnLinkLines(idx_CentroidXY_i,:);
        XYPtsOnLinkLinesByNuclei{i,1} = tmp_XYPtsOnLinkLines;
    end
end