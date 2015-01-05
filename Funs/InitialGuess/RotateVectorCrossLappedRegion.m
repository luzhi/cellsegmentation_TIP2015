function [ CentroidBoundaryPtsPairSet, XYPtsOnLinkLines ] = RotateVectorCrossLappedRegion( CentroidBoundaryPtsPairSet, XYPtsOnLinkLines )
% Fill Linking Lines in the Overlapping Region for Each Cell by Vector
% Rotation.
%
%   Legend:
%       CentroidBoundaryPtsPairSet{i,1} -- (x,y) of Centroid.
%       CentroidBoundaryPtsPairSet{i,2} -- (x,y) of Boundary.
%       
%       (x,y) - (column, row)
%   Usage: ...

%%
size_Original_CBPtsPair = size(CentroidBoundaryPtsPairSet,1);
%%
%+------------------------------------+
%| Find Junction Points for Each Cell |
%+------------------------------------+
JunctionPointSet = zeros(1,3);
numJunctionPoints = 0;
for i = 1:size(CentroidBoundaryPtsPairSet,1)
    %+--------------------------------+
    %|    End Boundary Points Case    |
    %+--------------------------------+
    if i == 1
        %+--------------------------------------+
        %|    If previous & following point     |
        %| are ALL different the the middle one |
        %+--------------------------------------+
        %|                  SKIP!               |
        %+--------------------------------------+
        if ~isequal(CentroidBoundaryPtsPairSet{i,1}, CentroidBoundaryPtsPairSet{size(CentroidBoundaryPtsPairSet,1),1}) ...
                && ~isequal(CentroidBoundaryPtsPairSet{i,1}, CentroidBoundaryPtsPairSet{(i + 1), 1})
            continue;
        end
        
        %+---------------------------------+
        %|       Is Junction Point?        |
        %+---------------------------------+
        %|          Either same as         |
        %| previous or following, not both |
        %+---------------------------------+
        if ~isequal(CentroidBoundaryPtsPairSet{i,1}, CentroidBoundaryPtsPairSet{size(CentroidBoundaryPtsPairSet,1),1}) ...
          || ~isequal(CentroidBoundaryPtsPairSet{i,1}, CentroidBoundaryPtsPairSet{(i + 1), 1})
      
            numJunctionPoints = numJunctionPoints + 1;

            JunctionPointSet(numJunctionPoints, 1) = CentroidBoundaryPtsPairSet{i,1}(1,1);  % Centroid Point - x
            JunctionPointSet(numJunctionPoints, 2) = CentroidBoundaryPtsPairSet{i,1}(1,2);  % Centroid Point - y
            JunctionPointSet(numJunctionPoints, 3) = CentroidBoundaryPtsPairSet{i,2}(1,1);  % Boundary Point - x
            JunctionPointSet(numJunctionPoints, 4) = CentroidBoundaryPtsPairSet{i,2}(1,2);  % Boundary Point - y
            JunctionPointSet(numJunctionPoints, 5) = i;     % Index in Centroid-Boundary Points Pairs Set
        end
    end
    
    if i == size(CentroidBoundaryPtsPairSet,1)
        %+--------------------------------------+
        %|    If previous & following point     |
        %| are ALL different the the middle one |
        %+--------------------------------------+
        %|                  SKIP!               |
        %+--------------------------------------+
        if ~isequal(CentroidBoundaryPtsPairSet{i,1}, CentroidBoundaryPtsPairSet{(i - 1), 1}) ...
                && ~isequal(CentroidBoundaryPtsPairSet{i,1}, CentroidBoundaryPtsPairSet{1,1})
            continue;
        end
        
        %+---------------------------------+
        %|       Is Junction Point?        |
        %+---------------------------------+
        %|          Either same as         |
        %| previous or following, not both |
        %+---------------------------------+
        if ~isequal(CentroidBoundaryPtsPairSet{i,1}, CentroidBoundaryPtsPairSet{(i - 1), 1}) ...
          || ~isequal(CentroidBoundaryPtsPairSet{i,1}, CentroidBoundaryPtsPairSet{1,1})
          
            numJunctionPoints = numJunctionPoints + 1;

            JunctionPointSet(numJunctionPoints, 1) = CentroidBoundaryPtsPairSet{i,1}(1,1);  % Centroid Point - x
            JunctionPointSet(numJunctionPoints, 2) = CentroidBoundaryPtsPairSet{i,1}(1,2);  % Centroid Point - y
            JunctionPointSet(numJunctionPoints, 3) = CentroidBoundaryPtsPairSet{i,2}(1,1);  % Boundary Point - x
            JunctionPointSet(numJunctionPoints, 4) = CentroidBoundaryPtsPairSet{i,2}(1,2);  % Boundary Point - y
            JunctionPointSet(numJunctionPoints, 5) = i;     % Index in Centroid-Boundary Points Pairs Set
        end
    end
    
    %+--------------------------------+
    %|      Internal Points Case      |
    %+--------------------------------+
    if i > 1 && i < size(CentroidBoundaryPtsPairSet,1)
        %+--------------------------------------+
        %|    If previous & following point     |
        %| are ALL different the the middle one |
        %+--------------------------------------+
        %|                  SKIP!               |
        %+--------------------------------------+
        if ~isequal(CentroidBoundaryPtsPairSet{i,1}, CentroidBoundaryPtsPairSet{(i - 1), 1}) ...
                && ~isequal(CentroidBoundaryPtsPairSet{i,1}, CentroidBoundaryPtsPairSet{(i + 1), 1})
            continue;
        end
        
        %+---------------------------------+
        %|       Is Junction Point?        |
        %+---------------------------------+
        %|          Either same as         |
        %| previous or following, not both |
        %+---------------------------------+
        if ~isequal(CentroidBoundaryPtsPairSet{i,1}, CentroidBoundaryPtsPairSet{(i - 1), 1}) ...
          || ~isequal(CentroidBoundaryPtsPairSet{i,1}, CentroidBoundaryPtsPairSet{(i + 1), 1})
            numJunctionPoints = numJunctionPoints + 1;
            
            JunctionPointSet(numJunctionPoints, 1) = CentroidBoundaryPtsPairSet{i,1}(1,1);  % Centroid Point - x
            JunctionPointSet(numJunctionPoints, 2) = CentroidBoundaryPtsPairSet{i,1}(1,2);  % Centroid Point - y
            JunctionPointSet(numJunctionPoints, 3) = CentroidBoundaryPtsPairSet{i,2}(1,1);  % Boundary Point - x
            JunctionPointSet(numJunctionPoints, 4) = CentroidBoundaryPtsPairSet{i,2}(1,2);  % Boundary Point - y
            JunctionPointSet(numJunctionPoints, 5) = i;     % Index in Centroid-Boundary Points Pairs Set
        end
    end
end

%%
%+-----------------------------------------------+
%| Check if size(JunctionPointSet,1) is (0,0,0)? |
%+-----------------------------------------------+
if size(JunctionPointSet,1) == 1
    if isequal(JunctionPointSet(1,:), [0,0,0])
        fprintf('No Junction Point here!\n');
    end    
else
    %%
    %+-------------------------------------------+
    %| Check if size(JunctionPointSet,1) is ODD? |
    %+-------------------------------------------+
    if mod(size(JunctionPointSet,1), 2) ~= 0
        %+-----------------------+
        %| Remove the ODD points |
        %+-----------------------+
        i = 1;
        while(i <= size(JunctionPointSet,1))
            %+-------------------------------------------+
            %| Process the Last Group of Points are ODD. |
            %+-------------------------------------------+
            if i == size(JunctionPointSet,1)
                JunctionPointSet(i,:) = [];
                continue;
            end
            
            %+-------------------------------------------------------+
            %| Process the Groups of Points other than the Last One. |
            %+-------------------------------------------------------+

            id_i = JunctionPointSet(i,end);
            id_i_plus_1 = JunctionPointSet(i + 1, end);

            %+--------------------+
            %| i is the ODD point |
            %+--------------------+
            if id_i_plus_1 - id_i > 1
                JunctionPointSet(i,:) = [];
                i = 1;
                continue;
            end
            
            i = i + 2;
        end
    end
    %%
    %+--------------------------------+
    %| Find Vectors Pair in Each Cell |
    %+--------------------------------+
    assert(mod(size(JunctionPointSet,1), 2) == 0, ...
            'Number of Junction Points MUST be even!');
        
        

        %+-----------------------------------------

    JunctionPointsPairSet = zeros(1, 6); %zeros((size(JunctionPointSet,1) / 2), 6);
    numJunctionPtsPairs = 0;
    VisitedBoundaryPtsSet = zeros(1,1);
    numVisitedBoundaryPts = 0;

    %+-----------------------------------------+
    %| Number of Boundary Points Pair is EVEN. |
    %+-----------------------------------------+
    %|    (Find the pair of Boundary Points    |
    %|  with same Centroid, positive/negative  |
    %|                  sequence)              |
    %+-----------------------------------------+
    if mod((size(JunctionPointSet,1) / 2),2) == 0
        for i = 1:2:(size(JunctionPointSet,1) - 1)
            XYcentroid_i_1 = JunctionPointSet(i,1:2);   % Centroid (x,y) of i
            XYcentroid_i_2 = JunctionPointSet(i + 1, 1:2);  % Centroid (x,y) of i + 1

            for j = 1:2:(size(JunctionPointSet,1)-1)
                if i == j
                    continue;
                end

                XYcentroid_j_1 = JunctionPointSet(j,1:2);   % Centroid (x,y) of j
                XYcentroid_j_2 = JunctionPointSet(j + 1, 1:2);  % Centroid (x,y) of j + 1

                if isequal(XYcentroid_i_1, XYcentroid_j_2) && isequal(XYcentroid_i_2, XYcentroid_j_1)
                    %+------------------------------+
                    %|      Check if Visited?       |
                    %+------------------------------+
                    if ismember(j, VisitedBoundaryPtsSet)
                        continue;
                    end

                    %+-------------------------------+
                    %| Shipping Boundary Points Pair |
                    %+-------------------------------+
                    numJunctionPtsPairs = numJunctionPtsPairs + 1;
                    JunctionPointsPairSet(numJunctionPtsPairs,1:2) = XYcentroid_i_1;        % Centroid (x,y) of i and j + 1
                    JunctionPointsPairSet(numJunctionPtsPairs,3:4) = JunctionPointSet(i, 3:4);    % Boundary point of i
                    JunctionPointsPairSet(numJunctionPtsPairs,5:6) = JunctionPointSet(j + 1, 3:4);    % Bounday point of j + 1
                    JunctionPointsPairSet(numJunctionPtsPairs,7) = JunctionPointSet(i, 5);
                    JunctionPointsPairSet(numJunctionPtsPairs,8) = JunctionPointSet(j + 1, 5);
                    
                    numJunctionPtsPairs = numJunctionPtsPairs + 1;
                    JunctionPointsPairSet(numJunctionPtsPairs, 1:2) = XYcentroid_i_2;   % Centroid (x,y) of i + 1 and j
                    JunctionPointsPairSet(numJunctionPtsPairs, 3:4) = JunctionPointSet(i + 1, 3:4);  % Boundary point of i + 1
                    JunctionPointsPairSet(numJunctionPtsPairs, 5:6) = JunctionPointSet(j, 3:4);      % Boundary point of j
                    JunctionPointsPairSet(numJunctionPtsPairs,7) = JunctionPointSet(i + 1, 5);
                    JunctionPointsPairSet(numJunctionPtsPairs,8) = JunctionPointSet(j, 5);

                    %+-----------------------------------------+
                    %| Checking-list for Visitied Points Pairs |
                    %|      (Positive/Negative Sequence)       |
                    %+-----------------------------------------+
                    numVisitedBoundaryPts = numVisitedBoundaryPts + 1;
                    VisitedBoundaryPtsSet(numVisitedBoundaryPts,1) = i;
                    continue;
                end
            end    
        end
    end

    %+-----------------------------------------+
    %| Number of Boundary Points Pair is EVEN. |
    %|  (But JunctionPointsPairSet is EMPTY!)  |
    %+-----------------------------------------+
    %|    (Find the farthest Boundary point    |
    %|          with same Centroid.)           |
    %+-----------------------------------------+
    if mod((size(JunctionPointSet,1) / 2),2) == 0
        if (max(JunctionPointsPairSet(:)) == 0) || (size(JunctionPointsPairSet,1) ~= (size(JunctionPointsPairSet,1) / 2))

            JunctionPointsPairSet = zeros(1, 6); %zeros((size(JunctionPointSet,1) / 2), 6);
            numJunctionPtsPairs = 0;
            VisitedBoundaryPtsSet = zeros(1,1);
            numVisitedBoundaryPts = 0;

            for i = 1:size(JunctionPointSet,1)
                if ismember(i, VisitedBoundaryPtsSet)
                    continue;
                end

                XYcentroid_i = JunctionPointSet(i,1:2);   % Centroid (x,y) of i

                %+----------------------------+
                %| Farthest Boundary Point ID |
                %+----------------------------+
                idx_SameCentroid_i = ismember(JunctionPointSet(:, 1:2), XYcentroid_i, 'rows');
                idx_ONES = find(idx_SameCentroid_i == 1);
                idx_visited_i = ismember(idx_ONES, VisitedBoundaryPtsSet);        
                idx_SameCentroid_i(idx_ONES(idx_visited_i == 1)) = 0;

                id_DistantBoundaryPt = find(idx_SameCentroid_i == 1, 1, 'last');

                %+-----------------------------+
                %|   Visited Boundary Points   |
                %+-----------------------------+
                %|            Skip!            |
                %+-----------------------------+
                if isempty(id_DistantBoundaryPt)
                    continue;
                end

                %+-------------------------------+
                %| Shipping Boundary Points Pair |
                %+-------------------------------+
                numJunctionPtsPairs = numJunctionPtsPairs + 1;
                JunctionPointsPairSet(numJunctionPtsPairs,1:2) = XYcentroid_i;  % Centroid (x,y) of i
                JunctionPointsPairSet(numJunctionPtsPairs,3:4) = JunctionPointSet(i,3:4);   % Boundary Point of i
                JunctionPointsPairSet(numJunctionPtsPairs,5:6) = JunctionPointSet(id_DistantBoundaryPt, 3:4);   % Boundary Point of farthest Boundary Point
                JunctionPointsPairSet(numJunctionPtsPairs,7) = JunctionPointSet(i, 5);
                JunctionPointsPairSet(numJunctionPtsPairs,8) = JunctionPointSet(id_DistantBoundaryPt, 5);

                %+-----------------------------------------+
                %| Checking-list for Visitied Points Pairs |
                %+-----------------------------------------+
                numVisitedBoundaryPts = size(VisitedBoundaryPtsSet,1) + 1;
                VisitedBoundaryPtsSet(numVisitedBoundaryPts, 1) = i;
                VisitedBoundaryPtsSet(numVisitedBoundaryPts + 1, 1) = id_DistantBoundaryPt;
            end
        end
    end
    %+----------------------------------------+
    %| Number of Boundary Points Pair is ODD. |
    %+----------------------------------------+
    %|    (Find the farthest Boundary point   |
    %|          with same Centroid.)          |
    %+----------------------------------------+
    if mod((size(JunctionPointSet,1) / 2), 2) == 1
        for i = 1:size(JunctionPointSet,1)
            if ismember(VisitedBoundaryPtsSet, i)
                continue;
            end

            XYcentroid_i = JunctionPointSet(i,1:2);   % Centroid (x,y) of i

            %+----------------------------+
            %| Farthest Boundary Point ID |
            %+----------------------------+
            idx_SameCentroid_i = ismember(JunctionPointSet(:, 1:2), XYcentroid_i, 'rows');
            idx_ONES = find(idx_SameCentroid_i == 1);
            idx_visited_i = ismember(idx_ONES, VisitedBoundaryPtsSet);        
            idx_SameCentroid_i(idx_ONES(idx_visited_i == 1)) = 0;

            id_DistantBoundaryPt = find(idx_SameCentroid_i == 1, 1, 'last');

            %+-----------------------------+
            %|   Visited Boundary Points   |
            %+-----------------------------+
            %|            Skip!            |
            %+-----------------------------+
            if isempty(id_DistantBoundaryPt)
                continue;
            end

            %+-------------------------------+
            %| Shipping Boundary Points Pair |
            %+-------------------------------+
            numJunctionPtsPairs = numJunctionPtsPairs + 1;
            JunctionPointsPairSet(numJunctionPtsPairs,1:2) = XYcentroid_i;  % Centroid (x,y) of i
            JunctionPointsPairSet(numJunctionPtsPairs,3:4) = JunctionPointSet(i,3:4);   % Boundary Point of i
            JunctionPointsPairSet(numJunctionPtsPairs,5:6) = JunctionPointSet(id_DistantBoundaryPt, 3:4);   % Boundary Point of farthest Boundary Point
            JunctionPointsPairSet(numJunctionPtsPairs,7) = JunctionPointSet(i, 5);
            JunctionPointsPairSet(numJunctionPtsPairs,8) = JunctionPointSet(id_DistantBoundaryPt, 5);
            
            %+-----------------------------------------+
            %| Checking-list for Visitied Points Pairs |
            %+-----------------------------------------+
            numVisitedBoundaryPts = size(VisitedBoundaryPtsSet,1) + 1;
            VisitedBoundaryPtsSet(numVisitedBoundaryPts, 1) = i;
            VisitedBoundaryPtsSet(numVisitedBoundaryPts + 1, 1) = id_DistantBoundaryPt;
        end
    end

    %%
    %+-----------------------------+
    %|      Vector Rotation        |
    %+-----------------------------+-----------+
    %|Structure of "JunctionPointsPairSet":    |
    %|      (Col 1, Col 2) [Centroid],         |
    %|      (Col 3, Col 4) [Boundary Point 1], |
    %|      (Col 5, Col 6) [Bounary Point 2],  |
    %+-----------------------------------------+
    
    %+-------------------------------+
    %|      Get Centroids (x,y)      |
    %+-------------------------------+
    CentroidsXY = zeros(1,1);
    numCentroids = 1;
    for i = 1:size(JunctionPointsPairSet,1)
        if ismember(JunctionPointsPairSet(i,1:2), CentroidsXY)
            continue;
        end
        
        CentroidsXY(numCentroids, 1:2) = JunctionPointsPairSet(i,1:2);
        numCentroids = numCentroids + 1;
    end
    
    %+------------------------------------+
    %| Group Boundary Points by Centroids |
    %+------------------------------------+
    [GroupCBPts, ~] = GroupCBPtsPairsByCentroids(CentroidBoundaryPtsPairSet, XYPtsOnLinkLines, CentroidsXY);
    
    %+-----------------------------------+
    %|      Auxiliary Match Method       |
    %|        for Vector Rotation        |
    %+-----------------------------------+
    for i = 1:size(GroupCBPts,1)
        GroupCBPts_i = GroupCBPts{i,1};
        Centroidxy_i = GroupCBPts_i{1,1};
        
        %+-------------------------------------+
        %|     Check if All Boundary Points    |
        %| in GroupCBPts_i is the centroidxy_i |
        %+-------------------------------------+
        CB_Same_list = zeros(size(GroupCBPts_i,1),1);
        for j = 1:size(GroupCBPts_i,1)
            if isequal(GroupCBPts_i{j,2}, Centroidxy_i)
                CB_Same_list(j,1) = 1;
            end
        end
        
        %+------------------------------------------+
        %|           Don't Rotate Vectors           |
        %| if all Boundary Points are the Centroids |
        %+------------------------------------------+
        idx_SamePts = find(CB_Same_list == 1);
        if size(idx_SamePts,1) == size(GroupCBPts_i,1)
            fprintf('Cell %d has crashed into the nuclei!\n', i);
            continue;
        end
        
        %+----------------------------+
        %| Junction Points for Cell i |
        %+----------------------------+
        jckPt_i = zeros(1,2);
        numJckPt = 1;
        for j = 1:size(JunctionPointsPairSet,1)
            if isequal(Centroidxy_i, JunctionPointsPairSet(j,1:2))
                jckPt_i(numJckPt,1:2) = JunctionPointsPairSet(j, 3:4);
                jckPt_i(numJckPt + 1, 1:2) = JunctionPointsPairSet(j, 5:6);
                numJckPt = numJckPt + 2;
            end
        end
        
        %+----------------------------+
        %|      Auxiliary Array       |
        %+----------------------------+
        Aux = -1 * ones(360, 3);
        Aux(:, 2:3) = 0;
        
        %+----------------------------+
        %|     Horizontal iVector     |
        %+----------------------------+
        iVec_hori = [1,0];
        
        for j = 1:size(GroupCBPts_i,1)         
            iVec_j = [GroupCBPts_i{j,2}(1,1) - Centroidxy_i(1,1), ...
                      GroupCBPts_i{j,2}(1,2) - Centroidxy_i(1,2)];
                  
            if isequal(iVec_j, [0,0])
                %+-------------------------+
                %|       Zero Vector       |
                %+-------------------------+
                if ismember(GroupCBPts_i{j,2}, jckPt_i, 'rows')
                    theta_j = 1;
                else
                    continue;
                end
            else
                %+-------------------------+
                %|     Non-zero Vector     |
                %+-------------------------+
                iVec_j = iVec_j / sqrt(sum(iVec_j .^ 2,2));

                cosine_hori_ivec_j = dot(iVec_hori, iVec_j);
                theta_j = acos(cosine_hori_ivec_j);

                %+---------------------------------------+
                %| Verify the theta_j or (2Pi - theta_j) |
                %+---------------------------------------+
                x_j = iVec_j(1,1);
                y_j = iVec_j(1,2);

                if x_j >= 0 && y_j >= 0
                    theta_j = round(theta_j * 180 / pi);
                elseif x_j < 0 && y_j >=0
                    theta_j = round(theta_j * 180 / pi);
                elseif x_j < 0 && y_j < 0
                    theta_j = round((2 * pi - theta_j) * 180 / pi);
                elseif x_j >= 0 && y_j < 0
                    theta_j = round((2 * pi - theta_j) * 180 / pi);
                end

                if theta_j == 360
                    theta_j = 359;
                end
            end            

            %+------------------------------+
            %|      Put iVec_j to Aux       |
            %+------------------------------+
            if ismember(Aux(theta_j + 1, 2:3), jckPt_i, 'rows')
                %//////////
                if ismember(GroupCBPts_i{j,2}, jckPt_i, 'rows')
                    isPut = 0;
                    k = theta_j;
                    while(isPut == 0)
                        if k < 1
                            k = 360;
                        end

                        if ~ismember(Aux(k, 2:3), jckPt_i, 'rows')
                            Aux(k,1) = 1;
                            Aux(k, 2:3) = GroupCBPts_i{j,2};
                            isPut = 1;
                        else
                            k = k - 1;
                        end
                    end
                else
                    continue;
                end
            else            
                Aux(theta_j + 1, 1) = 1;
                Aux(theta_j + 1, 2:3) = GroupCBPts_i{j,2};
            end
        end
        
        %+-------------------------------------------------------+
        %| Find angle between pair of Junction Points for Cell i |
        %+-------------------------------------------------------+
        for j = 1:size(JunctionPointsPairSet,1)
            if isequal(Centroidxy_i, JunctionPointsPairSet(j,1:2))
                jckPt_i_1 = JunctionPointsPairSet(j, 3:4);
                jckPt_i_2 = JunctionPointsPairSet(j, 5:6);
                
                memIndex_1 = ismember(Aux(:,2:3), jckPt_i_1, 'rows');
                memIndex_2 = ismember(Aux(:,2:3), jckPt_i_2, 'rows');
                
                id_jckPt_1 = find(memIndex_1 == 1);
                id_jckPt_2 = find(memIndex_2 == 1);
                
                id_start = min(id_jckPt_1, id_jckPt_2);
                id_end = max(id_jckPt_1, id_jckPt_2);
                
                num_of_one_1 = 0;
                for k = id_start:id_end
                    if Aux(k,1) == 1
                        num_of_one_1 = num_of_one_1 + 1;
                    end
                end
                angle1 = id_end - id_start;
                
                num_of_one_2 = 0;
                isStop = 0;
                id_aux = id_end;
                angle2 = 0;
                while(isStop == 0)
                    if id_aux >= 360
                        id_aux = 1;
                    end
                    
                    if id_aux == id_start;
                        isStop = 1;
                    end                    
                    
                    if Aux(id_aux,1) == 1
                        num_of_one_2 = num_of_one_2 + 1;
                    end
                    
                    id_aux = id_aux + 1;
                    angle2 = angle2 + 1;
                end
                
                if num_of_one_1 <= num_of_one_2
                    angle = angle1;
                else
                    angle = angle2;
                end
                
                %+------------------------------+
                %|  Find Start Vector by angle  |
                %+------------------------------+
                M = [cos(angle * pi / 180), -sin(angle * pi / 180); ...
                     sin(angle * pi / 180), cos(angle * pi / 180)];
                
                iVec_1 = [jckPt_i_1(1,1) - Centroidxy_i(1,1); ...
                          jckPt_i_1(1,2) - Centroidxy_i(1,2)];
                mag_1 = sqrt(sum(iVec_1 .^ 2, 1));
                iVec_1 = iVec_1 / mag_1;
                
                iVec_2 = [jckPt_i_2(1,1) - Centroidxy_i(1,1); ...
                          jckPt_i_2(1,2) - Centroidxy_i(1,2)];
                mag_2 = sqrt(sum(iVec_2 .^ 2, 1));
                iVec_2 = iVec_2 / mag_2;
                
                if dot(M * iVec_1, iVec_2) >= dot(M * iVec_2, iVec_1)
                    iVec_start = iVec_1;
                    VectorMagnitudeStart = mag_1;
                    VectorMagnitudeEnd = mag_2;
                else
                    iVec_start = iVec_2;
                    VectorMagnitudeStart = mag_2;
                    VectorMagnitudeEnd = mag_1;
                end
                
                %+------------------------------+
                %|    Rotate Vector by angle    |
                %+------------------------------+
                if angle > 1
                    numLines = angle;
                else
                    numLines = 1;
                end

                rotated_CentroidBoundaryPts_i = cell(1,2);
                for k = 1:numLines
                    Matrix_RotateStep = [cos(k * (pi / 180) *1), -sin(k * (pi / 180) *1); ...
                                         sin(k * (pi / 180) *1),  cos(k * (pi / 180) *1)];
                    iVec_k = Matrix_RotateStep * iVec_start;

                    if VectorMagnitudeStart < VectorMagnitudeEnd
                        VectorMagnitude_k = VectorMagnitudeStart + ...
                                                  k * (VectorMagnitudeEnd - VectorMagnitudeStart) / numLines;
                    else
                        if VectorMagnitudeStart > VectorMagnitudeEnd
                            VectorMagnitude_k = VectorMagnitudeStart - ...
                                                      k * (VectorMagnitudeStart - VectorMagnitudeEnd) / numLines;
                        else
                            if VectorMagnitudeStart == VectorMagnitudeEnd
                                if k < floor(numLines / 2)
                                    VectorMagnitude_k = VectorMagnitudeStart + ...
                                                  k * (VectorMagnitudeEnd - VectorMagnitudeStart) / numLines;
                                else
                                    VectorMagnitude_k = VectorMagnitudeStart - ...
                                                      k * (VectorMagnitudeStart - VectorMagnitudeEnd) / numLines;
                                end
                            end
                        end
                    end

                    Vec_k = (VectorMagnitude_k .* iVec_k)';
                    rotated_CentroidBoundaryPts_i{k,1} = Centroidxy_i;      % Centroid x,y
                    rotated_CentroidBoundaryPts_i{k,2} = [Centroidxy_i(1,1) + Vec_k(1,1), ...      % Boundary Point x
                                                          Centroidxy_i(1,2) + Vec_k(1,2)];         % Boundary Point y
                end
                CentroidBoundaryPtsPairSet = [CentroidBoundaryPtsPairSet; rotated_CentroidBoundaryPts_i];
            end
        end
    end

    %%
    %+---------------------------------------+
    %|    Generate New "XYPtsOnLinkLines"    |
    %+---------------------------------------+
    Rotated_XYPtsOnLinkLines = cell(1,1);
    numRotatedXYPtsOnLinkLines = 0;
    for i = (size_Original_CBPtsPair + 1):size(CentroidBoundaryPtsPairSet,1)
        d_i = sqrt((CentroidBoundaryPtsPairSet{i,1}(1,1) - CentroidBoundaryPtsPairSet{i,2}(1,1))^2 + ...
                   (CentroidBoundaryPtsPairSet{i,1}(1,2) - CentroidBoundaryPtsPairSet{i,2}(1,2))^2);
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
        
        numRotatedXYPtsOnLinkLines = numRotatedXYPtsOnLinkLines + 1;
        Rotated_XYPtsOnLinkLines{numRotatedXYPtsOnLinkLines,1}(:,1) = x_PtsLink;
        Rotated_XYPtsOnLinkLines{numRotatedXYPtsOnLinkLines,1}(:,2) = y_PtsLink;
    end
    
    XYPtsOnLinkLines = [XYPtsOnLinkLines; Rotated_XYPtsOnLinkLines];
end
end
