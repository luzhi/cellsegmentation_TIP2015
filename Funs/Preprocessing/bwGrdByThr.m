function [ imBwGrd ] = bwGrdByThr( imGrd, thresh )
% [ imBwGrd ] = bwGrdByThr( imGrd, thresh )
%   Generate binary image of the morphological gradient map by a threshold.

    imBwGrd = imGrd;
    imBwGrd(imBwGrd >= thresh) = 255;
    imBwGrd(imBwGrd < thresh) = 0;
%     imBwGrd(imBwGrd ~= 0) = 255;
end