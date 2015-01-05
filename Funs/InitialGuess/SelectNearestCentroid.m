function [CentroidBoundaryPtsPairSet, InternalCentroidIndicator, MinimalPartitionValue, XYPtsOnLinkLines] = ...
                SelectNearestCentroid(CentroidSet, ClumpBoundaryPts, ClumpMask)
%

% The (x,y) in this function is equivalent to that of mathematics (columns, rows).


    CentroidBoundaryPtsPairSet = cell(size(ClumpBoundaryPts,1),2);

    %%
    DistanceSet = cell(size(ClumpBoundaryPts,1),1);
    XY_PtsLinkSet = cell(size(ClumpBoundaryPts,1),1);
    %+------------------------------------------+
    %| Compute Distances & Pts on Linking Lines |
    %+------------------------------------------+
    for i = 1:size(ClumpBoundaryPts,1)    
        DistanceSet{i,1} = zeros(size(CentroidSet,1),1);
        XY_PtsLinkSet{i,1} = cell(size(CentroidSet,1),1);

        %+--------------------------------+
        %|   Visit Every Boundary Pixel   |
        %+--------------------------------+
        xBoundary = ClumpBoundaryPts(i,1);
        yBoundary = ClumpBoundaryPts(i,2);    

        for j = 1:size(CentroidSet,1)
            %+------------------------------+
            %|   Visit Every Nuclei Pixel   |
            %+------------------------------+
            xCentroid = CentroidSet(j,1);
            yCentroid = CentroidSet(j,2);

            %+------------------------------+
            %|      Distance between        |
            %| Boundary Pt i and Centroid j |
            %+------------------------------+
            d_ij = sqrt(double((xCentroid-xBoundary)^2 + (yCentroid-yBoundary)^2));
            d_ij = round(d_ij);

            %+--------------------------------+
            %|   Points on the linking line   |
            %+--------------------------------+
            x_PtsLink = zeros(d_ij + 1,1);
            y_PtsLink = zeros(d_ij + 1,1);

            x_PtsLink(1,1) = xCentroid;
            y_PtsLink(1,1) = yCentroid;
            x_PtsLink(d_ij+1,1) = xBoundary;
            y_PtsLink(d_ij+1,1) = yBoundary;
            for k = 2 : d_ij
                x_PtsLink(k,1) = (xCentroid * (d_ij + 1 - k) + xBoundary * ( k - 1)) / d_ij;
                y_PtsLink(k,1) = (yCentroid * (d_ij + 1 - k) + yBoundary * ( k - 1)) / d_ij;
            end

            %+----------------------------------+
            %| Put the Distance value of a pair |
            %|   into the Whole Distance Set    |
            %|      for further Selection       |
            %+----------------------------------+
            DistanceSet{i,1}(j,1) = d_ij;

            XY_PtsLinkSet{i,1}{j,1}(:,1) = x_PtsLink;
            XY_PtsLinkSet{i,1}{j,1}(:,2) = y_PtsLink;
        end
    end
    %%
    %+------------------------------------------+
    %| Mark Non-Linked Boundary Points AS Empty |
    %+------------------------------------------+
    for i = 1:size(ClumpBoundaryPts,1)
        id_minLinkLine = find(DistanceSet{i,1} == min(DistanceSet{i,1}(:)));
        %+--------------------------------+
        %| One Shortest Linking Line Case |
        %+--------------------------------+
        min_LinkXYPts = XY_PtsLinkSet{i,1}{id_minLinkLine(1,1)};
        LinkingLineMask = false(size(ClumpMask));
        for j = 1:size(min_LinkXYPts,1)
            LinkingLineMask(round(min_LinkXYPts(j,2)), round(min_LinkXYPts(j,1))) = 1;
        end

        CONDITION = ClumpMask .* LinkingLineMask;
        numRemainLinePts = size(find(CONDITION == 1),1);
        %+------------------------+
        %|     IS inside Clump    |
        %+------------------------+
        if (numRemainLinePts / size(min_LinkXYPts,1)) >= 0.95
            CentroidBoundaryPtsPairSet{i,1} = [min_LinkXYPts(1,1), min_LinkXYPts(1,2)];     % Centroid Points
            CentroidBoundaryPtsPairSet{i,2} = [min_LinkXYPts(end,1), min_LinkXYPts(end,2)]; % Boundary Points

            continue;
        else
            %+------------------------------+
            %| IS intersect with Background |
            %+------------------------------+
            CentroidBoundaryPtsPairSet{i,1} = [];
            CentroidBoundaryPtsPairSet{i,2} = [min_LinkXYPts(end,1), min_LinkXYPts(end,2)];
        end
    end
    %%
    %+------------------------------------------+
    %| Process Shortest Linking Lines intersect |
    %|          with Background Case            |
    %+------------------------------------------+
    for i = 1:size(CentroidBoundaryPtsPairSet,1)
        %+----------------------------+
        %| Skip Non-Empty Points Pair |
        %+----------------------------+
        if ~isempty(CentroidBoundaryPtsPairSet{i,1})
            continue;
        end

        %+-----------------------------+
        %|  Process Empty Points Pair  |
        %+-----------------------------+
        if isempty(CentroidBoundaryPtsPairSet{i,1})
            %+------------------------+
            %|  Internal Points Pair  |
            %+------------------------+
            if i > 1 && i < size(CentroidBoundaryPtsPairSet,1)            

                CB_PtsPair_before_i = [];
                CB_PtsPair_after_i = [];

                %+---------------------------+
                %| Find non-Empty Neighbours |
                %+---------------------------+
                j = 1;
                Found_nonEmpty = 0;
                numPtsToBeforeCBPair = 0;
                while( Found_nonEmpty == 0 && (i - j) >= 1)
                    if ~isempty(CentroidBoundaryPtsPairSet{i - j, 1})
                        CB_PtsPair_before_i = CentroidBoundaryPtsPairSet{i - j, 1};
                        Found_nonEmpty = 1;
                    end
                    numPtsToBeforeCBPair = numPtsToBeforeCBPair + 1;
                    j = j + 1;
                end

                j = 1;
                Found_nonEmpty = 0;
                numPtsToAfterCBPair = 0;
                while( Found_nonEmpty == 0 && (i + j) <= size(CentroidBoundaryPtsPairSet,1))
                    if ~isempty(CentroidBoundaryPtsPairSet{i + j, 1})
                        CB_PtsPair_after_i = CentroidBoundaryPtsPairSet{i + j, 1};
                        Found_nonEmpty = 1;
                    end
                    numPtsToAfterCBPair = numPtsToAfterCBPair + 1;
                    j = j + 1;
                end            

% %                 assert(~isempty(CB_PtsPair_before_i) && ~isempty(CB_PtsPair_after_i), ...
% %                     '[SelectNearestCentroid.m] This is NOT cell image!');
                
                %+-------------------------------+
                %| NO non-Empty Neighbour before |
                %+-------------------------------+
                if isempty(CB_PtsPair_before_i)
                    CentroidBoundaryPtsPairSet{i,1} = CB_PtsPair_after_i;
%                     fprintf('[SelectNearestCentroid.m] Should be NOT Happen! NO non-Empty Neighbour before.\n')
                end

                %+------------------------------+
                %| NO non-Empty Neighbour after |
                %+------------------------------+
                if isempty(CB_PtsPair_after_i)
                    CentroidBoundaryPtsPairSet{i,1} = CB_PtsPair_before_i;
                    continue;
%                     fprintf('[SelectNearestCentroid.m] Should be NOT Happen! NO non-Empty Neighbour after.\n')
                end

                %+------------------------------------------+
                %| Previous & Following Boundary Neighbours |
                %|          Has the Same Centroid,          |
                %|       The Empty Points Pair Links        |
                %|          to the Same Centroid            |
                %+------------------------------------------+
                if isequal(CB_PtsPair_before_i, CB_PtsPair_after_i)
                    CentroidBoundaryPtsPairSet{i,1} = CB_PtsPair_before_i;

                    continue;
                end

                %+------------------------------------------+
                %| Previous & Following Boundary Neighbours |
                %|         Has different Centroids          |
                %+------------------------------------------+
% %                 if ~isequal(CB_PtsPair_before_i, CB_PtsPair_after_i)
% %                     
% %                 end
                if ~isequal(CB_PtsPair_before_i, CB_PtsPair_after_i)
%                     d_to_before_i = sqrt((CentroidBoundaryPtsPairSet{i,2}(1,1) - CB_PtsPair_before_i(1,1))^2 + ...
%                                          (CentroidBoundaryPtsPairSet{i,2}(1,2) - CB_PtsPair_before_i(1,2))^2);
%                     d_to_after_i = sqrt((CentroidBoundaryPtsPairSet{i,2}(1,1) - CB_PtsPair_after_i(1,1))^2 + ...
%                                          (CentroidBoundaryPtsPairSet{i,2}(1,2) - CB_PtsPair_after_i(1,2))^2);

                    %+---------------------------------------------------+
                    %| Empty Points Pair and the Closest Points Pairs on |
                    %|        the Boundary has the same Centroid         |
                    %+---------------------------------------------------+
                    d_to_before_i = numPtsToBeforeCBPair;
                    d_to_after_i = numPtsToAfterCBPair;
                    
                    if d_to_before_i <= d_to_after_i
                        CentroidBoundaryPtsPairSet{i,1} = CB_PtsPair_before_i;
                    else
                        CentroidBoundaryPtsPairSet{i,1} = CB_PtsPair_after_i;
                    end
                    
                    continue;
                end
            end

            %+-------------------+
            %|  End Points Pair  |
            %+-------------------+
            if i == 1
                %+---------------------------+
                %| Find non-Empty Neighbours |
                %+---------------------------+
                j = 1;
                Found_nonEmpty = 0;
                while( Found_nonEmpty == 0 && (i + j) <= size(CentroidBoundaryPtsPairSet,1))
                    if ~isempty(CentroidBoundaryPtsPairSet{i + j, 1})
                        CB_PtsPair_after_i = CentroidBoundaryPtsPairSet{i + j, 1};
                        Found_nonEmpty = 1;
                    end
                    j = j + 1;
                end

                if Found_nonEmpty == 0
                    id_minLinkLine = find(DistanceSet{1,1} == min(DistanceSet{1,1}(:)));
                    min_LinkXYPts = XY_PtsLinkSet{1,1}{id_minLinkLine(1,1)};
                    CentroidBoundaryPtsPairSet{1,1} = [min_LinkXYPts(1,1), min_LinkXYPts(1,2)];     % Centroid Points
                    CentroidBoundaryPtsPairSet{1,2} = [min_LinkXYPts(end,1), min_LinkXYPts(end,2)]; % Boundary Points
                    continue;
                else
                    %+-------------------------------+
                    %| First Empty Points Pair Links |
                    %|    to the Centroid of its     |
                    %| Following non-Empty Neighbour |
                    %+-------------------------------+
                    CentroidBoundaryPtsPairSet{i,1} = CB_PtsPair_after_i;

                    continue;
                end
            end

            if i == size(CentroidBoundaryPtsPairSet,1)
                %+---------------------------+
                %| Find non-Empty Neighbours |
                %+---------------------------+
                j = 1;
                Found_nonEmpty = 0;
                while( Found_nonEmpty == 0 && (i - j) >= 1)
                    if ~isempty(CentroidBoundaryPtsPairSet{i - j, 1})
                        CB_PtsPair_before_i = CentroidBoundaryPtsPairSet{i - j, 1};
                        Found_nonEmpty = 1;
                    end
                    j = j + 1;
                end
                if Found_nonEmpty == 0
                    id_minLinkLine = find(DistanceSet{i,1} == min(DistanceSet{1,1}(:)));
                    min_LinkXYPts = XY_PtsLinkSet{i,1}{id_minLinkLine(1,1)};
                    CentroidBoundaryPtsPairSet{i,1} = [min_LinkXYPts(1,1), min_LinkXYPts(1,2)];     % Centroid Points
                    CentroidBoundaryPtsPairSet{i,2} = [min_LinkXYPts(end,1), min_LinkXYPts(end,2)]; % Boundary Points
                    continue;
                else
                    %+------------------------------+
                    %| Last Empty Points Pair Links |
                    %|    to the Centroid of its    |
                    %| Previous non-Empty Neighbour |
                    %+------------------------------+
                    CentroidBoundaryPtsPairSet{i,1} = CB_PtsPair_before_i;

                    continue;
                end
            end
        end
    end
    %%
    %+------------------------------------------------------------+
    %| Internal Cell does NOT have centroid-boundary points pairs |
    %+------------------------------------------------------------+

    %+-------------------------------+
    %|       Linked Centroids        |
    %+-------------------------------+
    CentroidCloseBoundary = zeros(size(CentroidBoundaryPtsPairSet,1),2);
    for i = 1:size(CentroidBoundaryPtsPairSet,1)
        CentroidCloseBoundary(i,1) = CentroidBoundaryPtsPairSet{i,1}(1,1);
        CentroidCloseBoundary(i,2) = CentroidBoundaryPtsPairSet{i,1}(1,2);
    end

    %+---------------------------+
    %|  Find Internal Centroids  |
    %+---------------------------+
    InternalCentroidIndicator = zeros(size(CentroidSet,1),1);
    for i = 1:size(CentroidSet,1)
        if ismember(CentroidSet(i,:), CentroidCloseBoundary)
            InternalCentroidIndicator(i,1) = 0;     % NOT Internal Cell
        else
            InternalCentroidIndicator(i,1) = 1;     % IS Internal Cell
        end
    end
    %%
    %+--------------------------+
    %|   Unique Distance List   |
    %+--------------------------+
    UniqueDistanceSet = zeros(size(CentroidBoundaryPtsPairSet,1),1);
    for i = 1:size(CentroidBoundaryPtsPairSet,1)
        UniqueDistanceSet(i,1) = sqrt((CentroidBoundaryPtsPairSet{i,1}(1,1) - CentroidBoundaryPtsPairSet{i,2}(1,1))^2 + ...
                                      (CentroidBoundaryPtsPairSet{i,1}(1,2) - CentroidBoundaryPtsPairSet{i,2}(1,2))^2);
    end
   
    assert(size(UniqueDistanceSet,1) == size(CentroidBoundaryPtsPairSet,1), ...
           'Distance Values Set MUST EQUAL to Points Pairs Set!');
    %%
    MinimalPartitionValue = 512 * 512;
    %+---------------------------------+
    %|   Points on Each Linking Line   |
    %+---------------------------------+
    XYPtsOnLinkLines = cell(size(CentroidBoundaryPtsPairSet,1),1);
    for i = 1:size(CentroidBoundaryPtsPairSet,1)
        d_i = UniqueDistanceSet(i,1);
        d_i = round(d_i);

        %+--------------------------------+
        %|   Points on the linking line   |
        %+--------------------------------+
        x_PtsLink = zeros(d_i + 1,1);
        y_PtsLink = zeros(d_i + 1,1);

        x_PtsLink(1,1) = CentroidBoundaryPtsPairSet{i,1}(1,1);   % Centroid X
        y_PtsLink(1,1) = CentroidBoundaryPtsPairSet{i,1}(1,2);   % Centroid Y
        x_PtsLink(d_i+1,1) = CentroidBoundaryPtsPairSet{i,2}(1,1);  % Boundary X
        y_PtsLink(d_i+1,1) = CentroidBoundaryPtsPairSet{i,2}(1,2);  % Boundary Y
        for j = 2 : d_i
            x_PtsLink(j,1) = (CentroidBoundaryPtsPairSet{i,1}(1,1) * (d_i + 1 - j) + ...
                              CentroidBoundaryPtsPairSet{i,2}(1,1) * ( j - 1)) / d_i;
            y_PtsLink(j,1) = (CentroidBoundaryPtsPairSet{i,1}(1,2) * (d_i + 1 - j) + ...
                              CentroidBoundaryPtsPairSet{i,2}(1,2) * ( j - 1)) / d_i;
        end

        XYPtsOnLinkLines{i,1}(:,1) = x_PtsLink;
        XYPtsOnLinkLines{i,1}(:,2) = y_PtsLink;
        
        %+-------------------------+
        %| Minimal Partition Value |
        %+-------------------------+
        if d_i <= MinimalPartitionValue
            MinimalPartitionValue = d_i;
        end
    end
end