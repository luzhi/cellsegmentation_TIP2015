function [ varargout ] = TestGMMSceneSegmentation( gmm_model_foreground, gmm_model_background, imSet, ID_TestImg )
% [ gmm_post ] = TestData4NCB_fullimage( gmm_model_nuclei, gmm_model_others, imSet, index_im2selected )
%
%   Compute the posterior probabilities (nuclei, clump and background) of
%   each pixel over the whole image using the GMM model of the two guys we
%   got from the 'TrainGMM2NCB(...)'.
%
%   Input:
%       gmm_model_clump: GMM model for clump (mean, cov, prior, etc.)
%       gmm_model_others: GMM model for background (mean, cov, prior, etc.)
%
%   Output:
%       gmm_post: a cell type of posterior probabilites of the two guys
%                 of each image
%           gmm_post{i,1}(2,:) - posterior of clump
%           gmm_post{i,1}(3,:) - posterior of background
%
%       ================================================================
%       Targets: P(y = Nuclei|X), P(y = Clump|X) and P(y = Background|X)
%       ================================================================

    % Parameters
%     imSize = 1024; %512;
    MAPPING = getmapping(8,'riu2');
    
    imTestFeaturesVector = [];
    for i = 1:length(ID_TestImg)
        im = imSet{ID_TestImg(i),1};
        [imSizeRow, imSizeCol] = size(im);
        imGradient = sobelgradient(im);
        J = lbp(im, 1, 8, MAPPING);
        imLBP = zeros(size(im));
        imLBP(2:(imSizeRow - 1),2:(imSizeCol - 1)) = J;
        
        %+------------------------+
        %| Vectorize Image Matrix |
        %+------------------------+
        VecImGray = reshape(im, 1, size(im,1)*size(im,2));
        VecImGradient = reshape(imGradient, 1, size(imGradient,1)*size(imGradient,2));
        VecImLBP = reshape(imLBP, 1, size(imLBP,1)*size(imLBP,2));
        
        imTestFeaturesVector.data{i,1} = [ double( VecImGray ) / 255; ...
            double( VecImGradient) / double( max(VecImGradient) ); ...
            double( VecImLBP / 255 ) ];
    end
    
    gmm_post = cell(length(ID_TestImg),1);
    for i = 1:length(ID_TestImg)
        % Test with training data
        prob_XY_foreground = pdfgauss( imTestFeaturesVector.data{i,1}, gmm_model_foreground );
        prob_XY_background = pdfgauss( imTestFeaturesVector.data{i,1}, gmm_model_background );    

        prior_Y_foreground = gmm_model_foreground.Prior;
        prior_Y_background = gmm_model_background.Prior;
        
        post_Y_foreground = (prob_XY_foreground * prior_Y_foreground) ./ ( prob_XY_foreground * prior_Y_foreground + prob_XY_background * prior_Y_background );
        post_Y_background = (prob_XY_background * prior_Y_background) ./ ( prob_XY_foreground * prior_Y_foreground + prob_XY_background * prior_Y_background );

        gmm_post{i,1} = zeros(2, imSizeRow * imSizeCol);
        gmm_post{i,1}(1, :) = post_Y_foreground(1,:);
        gmm_post{i,1}(2, :) = post_Y_background(1,:);
    end
    
    varargout{1,1} = gmm_post;
end

