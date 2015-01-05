function [ NucleiClumpAreaRatio ] = ComputeAreaRatio_NucleiClump(nucleiRegionAreaSize, nucleiRegionPixelIdxList, clumpMask)
%
%
%   Compute a set of area ratios between nuclei regions and its clump
%   regions
    
    NucleiClumpAreaRatio = -1;
    clumpStats = regionprops(clumpMask, 'Area', 'PixelIdxList');


    for k = 1:length(clumpStats)
        isNucleiBelongToClump = isequal(intersect(nucleiRegionPixelIdxList(:), clumpStats(k,1).PixelIdxList(:)), ...
                                        nucleiRegionPixelIdxList(:));
        if isNucleiBelongToClump == 1
            NucleiClumpAreaRatio = double( nucleiRegionAreaSize / clumpStats(k,1).Area );
            break;
        end
    end

    if NucleiClumpAreaRatio == -1;
        NucleiClumpAreaRatio = 0.000001;
    end
end