function [ varargout ] = TrainGMMSceneSegmentation( imSet, ForeBackgroundMaskSet, idx_TrainImg )
% [ gmm_model_nuclei, gmm_model_clump, gmm_model_background ] = TrainGMMForeBackground_noNuclei4Training( imSet, imGTSet, imClumpEdgesSet, index_im2train )
%   Train the parameters (MuK, SigmaK, PiK) of a GMM for the clump
%   (including nuclei) and the background (named as 'others')
%   
%   Features: gray + gradient + lbp
%
%   Input:
%       imSet: a set of images for training
%       imNCBMasks: a cell type of the masks of nuclei, clump and
%                   background
%                   nuclei: 0
%                   clump:  100
%                   background: 255
%       index_im2train: radomly selected images for training
%
%   Output:
%       gmm_model_clump: GMM model for clump (mean, cov, prior, etc.)
%       gmm_model_others: GMM model for background (mean, cov, prior, etc.)
%
%       ===============================================
%       Training targets: P(X|theta_C) and P(X|theta_B)
%       ===============================================
%
%   Example:
%       imSet = cell(4,1);
%       imSet{1,1} = imread('ims\EDF000.png');
%       imSet{2,1} = imread('ims\EDF001.png');
%       imSet{3,1} = imread('ims\EDF002.png');
%       imSet{4,1} = imread('ims\EDF003.png');
%       imGTSet = cell(4,1);
%       imGTSet{1,1} = imread('ims\EDF000_GTMask.png');
%       imGTSet{2,1} = imread('ims\EDF001_GTMask.png');
%       imGTSet{3,1} = imread('ims\EDF002_GTMask.png');
%       imGTSet{4,1} = imread('ims\EDF003_GTMask.png');
%
%       [ gmm_model_clump, gmm_model_others ] = ...
%       TrainGMMForeBackground_noNuclei4Training( imSet, imGTSet, imClumpEdgesSet, index_im2train );

    % Parameters
%     imSize = 1024; %512;
    MAPPING = getmapping(8,'riu2');

    tmpData4train = [];
    
    tmpData4train.foreground_gray = -1;
    tmpData4train.foreground_gradient = -1;
    tmpData4train.foreground_lbp = -1;
    
    tmpData4train.background_gray = -1;
    tmpData4train.background_gradient = -1;
    tmpData4train.background_lbp = -1;
    
    data4train = [];    
    
    for i = 1:size(idx_TrainImg,1)
        % features extraction
        im = imSet{idx_TrainImg(i),1};
        [imSizeRow, imSizeCol] = size(im);
        imForeBackgroundMask = ForeBackgroundMaskSet{idx_TrainImg(i),1};
        imGradient = sobelgradient(im);        
        J = lbp(im, 1, 8, MAPPING);
        imLBP = zeros(size(im));
        imLBP(2:(imSizeRow - 1),2:(imSizeCol - 1)) = J;
        
        % ROI of the image    
        
        idx_foreground = imForeBackgroundMask == 100;
        
        idx_background = imForeBackgroundMask == 255;
        
        tmpData4train.foreground_gray =[tmpData4train.foreground_gray; im( idx_foreground )];
        tmpData4train.foreground_gradient = [tmpData4train.foreground_gradient; imGradient( idx_foreground )];
        tmpData4train.foreground_lbp = [tmpData4train.foreground_lbp; imLBP( idx_foreground )];
        
        tmpData4train.background_gray = [tmpData4train.background_gray; im( idx_background )];
        tmpData4train.background_gradient = [tmpData4train.background_gradient; imGradient( idx_background )];
        tmpData4train.background_lbp = [tmpData4train.background_lbp; imLBP( idx_background )];
    end
    
    % the Training data set! + normalization
    
    data4train.foreground = [tmpData4train.foreground_gray(2:length(tmpData4train.foreground_gray)), ...
                        tmpData4train.foreground_gradient(2:length(tmpData4train.foreground_gradient)), ...
                        tmpData4train.foreground_lbp(2:length(tmpData4train.foreground_lbp))];
    data4train.foreground = double(data4train.foreground);
    data4train.foreground(:,1) = double( data4train.foreground(:,1) / 255 );
    data4train.foreground(:,2) = double( data4train.foreground(:,2) / max(data4train.foreground(:,2)) );
    data4train.foreground(:,3) = double( data4train.foreground(:,3) / 255 );
    
    data4train.foreground_labels = ones(length(tmpData4train.foreground_gray) - 1,1);
    
    data4train.background = [tmpData4train.background_gray(2:length(tmpData4train.background_gray)), ...
                         tmpData4train.background_gradient(2:length(tmpData4train.background_gradient)), ...
                         tmpData4train.background_lbp(2:length(tmpData4train.background_lbp))];
    data4train.background = double(data4train.background);
    data4train.background(:,1) = double( data4train.background(:,1) / 255 );
    data4train.background(:,2) = double( data4train.background(:,2) / max(data4train.background(:,2)) );
    data4train.background(:,3) = double( data4train.background(:,3) / 255 );
    
    data4train.background_labels = 2 * ones(length(tmpData4train.background_gray) - 1,1);

    DSTrain_foreground = [];
    DSTrain_foreground.X = data4train.foreground(:,:)';    
    DSTrain_foreground.y = data4train.foreground_labels(:,:)'; 
    DSTrain_foreground.name = 'Training set - Foreground';
    DSTrain_foreground.dim = size(data4train.foreground,2);
    DSTrain_foreground.num_data = size(data4train.foreground,1);
    
    DSTrain_background = [];
    DSTrain_background.X = data4train.background(:,:)';
    DSTrain_background.y = data4train.background_labels(:,:)';
    DSTrain_background.name = 'Training set - Background';
    DSTrain_background.dim = size(data4train.background,2);
    DSTrain_background.num_data = size(data4train.background,1);
    
    % Estimate the paramters of GMM
    gmm_model_foreground = mlcgmm(DSTrain_foreground);
    gmm_model_background = mlcgmm(DSTrain_background);
    
    % Ouput estimated GMM parameters
    varargout{1,1} = gmm_model_foreground;
    varargout{2,1} = gmm_model_background;
end

