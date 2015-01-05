function [ imGrd ] = morphGrad( im )
% [ imGrd ] = morphGrad( im )
%   Generate a morphological gradient map for an image


    SE = strel('disk',1,8);
    
    imGrd = imdilate(im,SE) - imerode(im,SE);
end

