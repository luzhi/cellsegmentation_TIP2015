function [ varargout ] = ComputeConfidenceAsScene( gmm_post, imRol, imCol )
% [ varargout ] = ComputeConfidenceNCB( gmm_post )
%   Depending on the posterior obtained by GMM, calculate the labels for
%   each pixels by confidence
%
%   Method: compute ratio between any pair of NCB probabilites
%
%   Input:
%       gmm_post: elements of a cell of posteriors of training images (NCB)
%
%   Output:
%       varargout{2,1} - confidence between clump and background


    im_foreground = reshape(gmm_post(1,:), imRol, imCol);
    im_background = reshape(gmm_post(2,:), imRol, imCol);
    
    % confidence between clump and background
    im_conf_foreground_background = exp( im_foreground - im_background );
    im_conf_foreground_background = im_conf_foreground_background / max(max( im_conf_foreground_background ));
    im_conf_foreground_background = uint8( im_conf_foreground_background * 255 );
    im_conf_foreground_background = im2bw( im_conf_foreground_background );
    
    varargout{1,1} = im_conf_foreground_background;
end

