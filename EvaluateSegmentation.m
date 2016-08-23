function [DSC, FNRo, TPRp, FDRo, stdDSC, stdFNRo, stdTPRp, stdFDRo] = ...
    EvaluateSegmentation(groundTruth, segmentationResult, diceRatio)
% EvaluateSegmentation calculates Dice Similarity Coefficient (or F-measure
% at pixel level), False Negative Rate (object level), True Positive Rate
% (pixel level), False Discovery Rate (object level) and their
% corresponding standard deviations.
% 
% Inputs:
%   1) groundTruth - ground truth segmentation of objects. It needs to be a
%   cell array (of size (M, 1), where M is the number of segmented images.
%   Also, each cell element is a cell array itself (of size (N_i, 1), where
%   N_i is the number of objects image i. groundTruth{i}{j} contains a mask
%   of segmentation of object j in image i.
% 
%   2) segmentationResult - segmentation result of objects. It needs to be
%   structured similar to groundTruch cell array. However, the number of
%   segmented objects in each image may be different.
% 
%   3) diceRatio - the minimum dice similarity coefficient of a "good"
%   segmentation. Default value is 0.7.
% 
% Outputs:
%   1) DSC - average dice similarity coefficient of all good segmentations
%   (in all images).
% 
%   2) FNRo - false negative rate at object level. It is equal to "the
%   total number of good segmentations" divided by "the total number of
%   objects in ground truth data".
% 
%   3) TPRp - average true positive rate (at pixel level) of all good
%   segmentations (in all images).
% 
%   4) FDRo - false discovery rate at object level. It is equal to "the
%   total number of segmented objects minus the total number of good
%   segmentations" divided by "the total number of objects in segmentation
%   result".
% 
%   5) stdDSC - standard deviation of average dice similarity coeffecients
%   of good segmentations in each image.
% 
%   6) stdFNRo - standard deviation of false negative rates at object level
%   of images.
% 
%   7) stdTPRp - standard deviation of true positive rates at pixel level 
%   of good segmentations in each image.
% 
%   8) stdFDRo - standard deviation of false discovery rates at object
%   level of images.
% 
% Example:
%	[DSC, FNRo, TPRp, FDRo, stdDSC, stdFNRo, stdTPRp, stdFDRo] = ...
%	    EvaluateSegmentation(groundTruth, segmentationResult);
%	fprintf(repmat(['\t%.3f', char(177), '%.3f'], 1, 4), ...
%	    DSC, stdDSC, FNRo, stdFNRo, TPRp, stdTPRp, FDRo, stdFDRo);
% 
% Copyright (c) 2016, Hady Ahmady Phoulady
% Department of Computer Science and Engineering,
% University of South Florida, Tampa, FL.
%
% Last modified: August 22, 2016

% Set default value for diceRatio in case it is not specificed
if (~exist('diceRatio', 'var'))
    diceRatio = 0.7;
end

[cellsDSC, cellsTPRp, allPairsDice, allPairsTPRp] = ...
    deal(cell(length(groundTruth), 1));
[imagesFNo, imagesFDo, averageImageDSC, imagesFNRo, imagesFDRo, ...
    averageImageTPRp] = deal(zeros(length(groundTruth), 1));

for i = 1: length(groundTruth)
    [allPairsDice{i}, allPairsTPRp{i}] = ...
        deal(zeros(length(groundTruth{i}), length(segmentationResult{i})));
    segArea = zeros(length(segmentationResult{i}), 1);

    % Iterate through images and calculate pairwise dice similarity
    % coefficient and true positive rate at pixel level of each object
    % in ground truth with each object in segmentation result
    for s = 1: length(segmentationResult{i})
        segArea(s) = nnz(segmentationResult{i}{s});
    end
    for g = 1: length(groundTruth{i})
        gtArea = nnz(groundTruth{i}{g});
        for s = 1: length(segmentationResult{i})
            intArea = nnz(groundTruth{i}{g} & segmentationResult{i}{s});
            allPairsDice{i}(g, s) = 2 * intArea / (gtArea + segArea(s));
            allPairsTPRp{i}(g, s) = intArea / gtArea;
        end
    end

    % Assign segmented regions to ground truth regions starting from
    % the highest pairwise dice similarity coefficient
    tmpAllPairsDice = allPairsDice{i};
    [allMaxDice, allMaxInd] = deal(zeros(length(groundTruth{i}), 1));
    origGTInd = (1: length(groundTruth{i}));
    origSInd = (1: length(segmentationResult{i}));
    for g = 1: min(length(groundTruth{i}), length(segmentationResult{i}))
        [maxDiceArr, maxIndArr] = max(tmpAllPairsDice, [], 2);
        [maxDice, maxGTInd] = max(maxDiceArr);
        maxSInd = maxIndArr(maxGTInd);
        allMaxDice(origGTInd(maxGTInd)) = maxDice;
        allMaxInd(origGTInd(maxGTInd)) = origSInd(maxSInd);
        tmpAllPairsDice(maxGTInd, :) = [];
        origGTInd(maxGTInd) = [];
        tmpAllPairsDice(:, maxSInd, :) = [];
        origSInd(maxSInd) = [];
    end

    % Compute the measures
    cellsDSC{i} = allMaxDice(allMaxDice > diceRatio);
    averageImageDSC(i) = mean(cellsDSC{i});
    cellsTPRp{i} = allPairsTPRp{i}(sub2ind(size(allPairsTPRp{i}), ...
        find(allMaxDice > diceRatio), allMaxInd(allMaxDice > diceRatio)));
    averageImageTPRp(i) = mean(cellsTPRp{i});
    imagesFNo(i) = nnz(allMaxDice <= diceRatio);
    imagesFNRo(i) = imagesFNo(i) / length(groundTruth{i});
    imagesFDo(i) = (length(segmentationResult{i}) - ...
        nnz(allMaxDice > diceRatio));
    imagesFDRo(i) = imagesFDo(i) / length(segmentationResult{i});
end

% Remove NaN elements that happen when all cells are missed, there is
% no reported cell or there is no cell in the image ground truth.
averageImageDSC(isnan(averageImageDSC)) = 0;
averageImageTPRp(isnan(averageImageTPRp)) = 0;
imagesFNRo(isnan(imagesFNRo)) = 0;
imagesFDRo(isnan(imagesFDRo)) = 0;

% Set the outputs
[DSC, FNRo, TPRp, FDRo, stdDSC, stdFNRo, stdTPRp, stdFDRo] = ...
    deal(mean(cell2mat(cellsDSC)), ...
    sum(imagesFNo) / sum(cellfun(@length, groundTruth)), ...
    mean(cell2mat(cellsTPRp)), ...
    sum(imagesFDo) / sum(cellfun(@length, segmentationResult)),  ...
    std(averageImageDSC),  ...
    std(imagesFNRo), std(averageImageTPRp), std(imagesFDRo));
end
    