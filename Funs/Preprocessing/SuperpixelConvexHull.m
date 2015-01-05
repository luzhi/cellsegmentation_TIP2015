function [ ConvexhullSceneMask ] = SuperpixelConvexHull( im,  spRATIO, spKERNELSIZE, spMAXDIST )
% Preprocessing - ...
% X,Y,K and V for segmentation over each convex hull
% K - indicator of the hulls
% V - the volume (number of pixels) in each hull
% tic;
    %%
    %+--------------------+
    %| Parameters Setting |
    %+--------------------+
    handles.ad_k = 0.05;
    handles.ad_iter = 5;
    handles.conhull_thr = 100;
    handles.thr_tinyfrags = 400;
    
    %%
    %+---------+
    %| Denoise |
    %+---------+
    img_denoised = anisodiff(im, handles.ad_iter, ...
                             handles.ad_k, 0.25, 1);
    img_denoised = uint8(img_denoised);
    
    %+---------------------+
    %| Compute Superpixels |
    %|    by QuickShift    |
    %+---------------------+
    [imSegbySP, junk] = vl_quickseg(img_denoised, spRATIO, spKERNELSIZE, spMAXDIST); %0.5,2,5);
    
    %+-------------------------+
    %| Morphorlogical Gradient |
    %+-------------------------+
    imGrdbyMorph = morphGrad(imSegbySP);
    bwImGrdbyMorph = bwGrdByThr(imGrdbyMorph,0.03);
    bwImGrdbyMorph = bwmorph(bwImGrdbyMorph,'thin');  
    
    %+----------------------------------------+
    %| Obtain Convex Hull for Each Superpixel |
    %+----------------------------------------+
    [SuperPixelLabelsMat, junk] = bwlabel(bwImGrdbyMorph, 8);
    % X,Y,K and V for segmentation over each convex hull
    % K - indicator of the hulls
    % V - the volume (number of pixels) in each hull

    ConvexhullSceneMask = false(size(im));
    for i=1:max(SuperPixelLabelsMat(:))
        [Y,X] = find(SuperPixelLabelsMat==i);
        if (length(Y) > handles.conhull_thr)            
            [K, junk] = convhull(X,Y);                        
            SuperpixelConvexhull_i = roipoly(im,X(K),Y(K));
            ConvexhullSceneMask( SuperpixelConvexhull_i == 1 ) = 1;
        end;
    end;
    
    ConvexhullSceneMask = im2bw(ConvexhullSceneMask);                
end

