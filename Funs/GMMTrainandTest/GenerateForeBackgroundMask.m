function [ ForeBackgroundMaskSet ] = GenerateForeBackgroundMask( ConvexhullSceneMaskSet )
% Assign 255 to the background region and 100 to the foreground region
% (i.e., nuclei or clump)
%   
%   Input:
%  --------
%     ConvexhullSceneMaskSet - Scene mask obtained by convex hull in Step 1
%   Output:
%  --------
%     ForeBackgroundMaskSet - Each element, 100 refers to scene object, and
%                             255 means background region.
    
    ForeBackgroundMaskSet = cell(size(ConvexhullSceneMaskSet,1),1);
    for i = 1:size(ConvexhullSceneMaskSet,1)
        
        %+-------------------------+
        %| Gray = 255 is the value |
        %|           of            |
        %|       Background        |
        %+-------------------------+
        imCanvas_i = 255 * ones(size(ConvexhullSceneMaskSet{i,1}));

        %+-------------------------+
        %| Gray = 100 is the value |
        %|           of            |
        %|       Foreground        |
        %+-------------------------+
        ForegroundMask_i = ConvexhullSceneMaskSet{i,1};
        imCanvas_i(sub2ind(size(imCanvas_i), find(ForegroundMask_i == 1))) = 100;

        ForeBackgroundMaskSet{i,1} = uint8(imCanvas_i);
    end
end

