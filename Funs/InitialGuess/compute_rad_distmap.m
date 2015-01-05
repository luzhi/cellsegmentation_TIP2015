function [InitialGuessSet] = compute_rad_distmap(ClumpMask, nucleusMaskInsideClump, beta)

[rowClumpMask,columnClumpMask]=size(ClumpMask);

%+-----------------------+
%| Clump Boundary Pixels |
%+-----------------------+
objBoundary = bwboundaries(ClumpMask);
BoundaryPixels(:,1) = objBoundary{1,1}(:,2);     % x-axis <--> column
BoundaryPixels(:,2) = objBoundary{1,1}(:,1);     % y-axis <--> row

%+----------------------------------------------+
%| Centroid-Boundary Points Pairs for ALL Cells |
%+----------------------------------------------+
%|       Boundary Cells or Internal Cells       |
%+----------------------------------------------+
[CentroidBoundaryPtsPairSet, MinimalPartitionValue, XYPtsOnLinkLines ] = ...
        ComputeCentroidBoundaryPairs( ClumpMask, nucleusMaskInsideClump );

% MinimalPartitionValue = round(MinimalPartitionValue);

%+-------------------------------------------------+
%| Group "CentroidBoundaryPtsPairSet" by Centroids |
%+-------------------------------------------------+---------------------+
%| "CentroidBoundaryPtsPairSetByNuclei" structure:                       |
%|      {[Centroid x, Centroid y]}{[Boundary_Point x, Boundary_Point y]} |
%+-----------------------------------------------------------------------+

%+-------------------+
%| Nucleus Centroids |
%+-------------------+
nucleusRegions = regionprops(nucleusMaskInsideClump, 'Centroid', 'PixelList');
CentroidsXY = zeros(length(nucleusRegions),2);
for i = 1:length(nucleusRegions)
    CentroidsXY(i,1) = nucleusRegions(i,1).Centroid(:,1);      % column - x
    CentroidsXY(i,2) = nucleusRegions(i,1).Centroid(:,2);      % row - y    
end

%+---------------------------------------+
%| Group "Centroid-Boundary Points Pairs |
%|            by Centroids               |
%+---------------------------------------+
[CentroidBoundaryPtsPairSetByNuclei, XYPtsOnLinkLinesByNuclei] = ...
    GroupCBPtsPairsByCentroids(CentroidBoundaryPtsPairSet, XYPtsOnLinkLines, CentroidsXY);

%+---------------------------------------------+
%| Assign Distance Values by Logistic Function |
%+---------------------------------------------+
%|  Draw the Initial Guess Mask for Each Cell  |
%+---------------------------------------------+

    InitialGuessSet = cell(size(CentroidsXY,1),2);  % Column 1: mask; Column 2: nuclei centroid (x,y)
    
    for i = 1:size(CentroidsXY,1)
        %+-------------------------------------------------+
        %| Compute distance values on Interpolation Points |
        %+-------------------------------------------------+
        CBPtsPairSetofaCell = CentroidBoundaryPtsPairSetByNuclei{i,1};
        radLinePartitionSet = cell(1,1);
        InitialGuessSet{i,1} = zeros(rowClumpMask,columnClumpMask);
        InitialGuessSet{i,2} = CBPtsPairSetofaCell{1,1};
        for j = 1 : size(CBPtsPairSetofaCell,1) 
            %+---------------------------+
            %| Interpolation Poins (x,y) |
            %|      for Each Line.       |
            %+---------------------------+
            InterpolationPointsXYOnLine = ...
                GetInterpolationPointsXY(CBPtsPairSetofaCell{j,1}, CBPtsPairSetofaCell{j,2}, MinimalPartitionValue);

            %+-----------------------------------------+
            %| Distance Values at Interpolation Points |
            %|              for Each Line.             |
            %+-----------------------------------------+
            %|      Logistic Function with "beta"      |
            %+-----------------------------------------+
            partitionDist = [0:1/MinimalPartitionValue:1]';
            partitionDist = -2 ./ (1 + exp(-beta .* partitionDist)) + 2;
            radLinePartitionSet{j,1} = [InterpolationPointsXYOnLine(:,1), InterpolationPointsXYOnLine(:,2), ...
                                        partitionDist(:,1)];
        end

        %+--------------------------------------+
        %| Interpolation for Each Initial Guess |
        %+--------------------------------------+
        interpretX = [];
        interpretY = [];
        interpretValue = [];
        for j = 1 : size(CBPtsPairSetofaCell,1)
            interpretY = [interpretY; radLinePartitionSet{j,1}(:,2)];
            interpretX = [interpretX; radLinePartitionSet{j,1}(:,1)];
            interpretValue = [interpretValue; radLinePartitionSet{j,1}(:,3)];
        end      
        
        if isempty(interpretX) || isempty(interpretY) || isempty(interpretValue)
            interpretY = [0];
            interpretX = [0];
            interpretValue = [0];
        end
      
        F = [];
        F = TriScatteredInterp(interpretX,interpretY,interpretValue);

        for j = 1 : rowClumpMask
            for k = 1 : columnClumpMask                
                if(isnan(F(k,j))) || (isempty(F(k,j)))
                    InitialGuessSet{i,1}(j,k) = 0;
                else
                    InitialGuessSet{i,1}(j,k) = F(k,j);
                end
            end
        end
        
        %+-----------------------------+
        %| Remove values outside clump |
        %+-----------------------------+
        InitialGuessSet{i,1} = InitialGuessSet{i,1} .* ClumpMask;        
       
        
        %+-------------------------------+
        %| Force values at boundary as 0 |
        %+-------------------------------+
        for j = 1:size(BoundaryPixels,1)
            InitialGuessSet{i,1}(BoundaryPixels(j,2), BoundaryPixels(j,1)) = 0;
        end
        
        %+------------------------------------+
        %| Assign lagest values outside clump |
        %+------------------------------------+
        InitialGuessSet{i,1}(ClumpMask ~= 1) = 1;
    end
end