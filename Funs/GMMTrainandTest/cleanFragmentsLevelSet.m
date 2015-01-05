function [ varargout ] = cleanFragmentsLevelSet( im, ForeBackgroundMASK, iter_inner, iter_outer, alfa, lambda )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
    
    Img=double(im(:,:,1));
    
    %+-------------------------+
    %| Fixed Parameter Setting |
    %+-------------------------+
    timestep=5;  % time step
    mu=0.2/timestep;  % coefficient of the distance regularization term R(phi)
    epsilon=1.5; % papramater that specifies the width of the DiracDelta function

    %+----------------------+
    %| Compute Edge Map "g" |
    %+----------------------+
    sigma=1.5;     % scale parameter in Gaussian kernel
    G=fspecial('gaussian',15,sigma);
    Img_smooth=conv2(Img,G,'same');  % smooth image by Gaussiin convolution
    [Ix,Iy]=gradient(Img_smooth);
    f=Ix.^2+Iy.^2;
    g=1./(1+f);  % edge indicator function.

    %+----------------------------------------+
    %| Initialize LSF as binary step function |
    %|      Object = -c0, Background = c0     |
    %+----------------------------------------+
    c0 = 2;
    initialLSF = c0*ones(size(Img));
    initialLSF(ForeBackgroundMASK == 1) = -c0;  
    phi = initialLSF;

    potential=2;  
    if potential == 1
        potentialFunction = 'single-well';  % use single well potential p1(s)=0.5*(s-1)^2, which is good for region-based model 
    elseif potential == 2
        potentialFunction = 'double-well';  % use double-well potential in Eq. (16), which is good for both edge and region based models
    else
        potentialFunction = 'double-well';  % default choice of potential function
    end


    %+---------------------------+
    %| Start Level Set Evolution |
    %+---------------------------+
    for n=1:iter_outer
        fprintf('\tRemove Fragments by Level Set - Iteration #: %d\n', n);
        phi = drlse_Denoise(phi, g, lambda, mu, alfa, epsilon, timestep, iter_inner, potentialFunction);
    end

    %+----------------------------+
    %| Refine Level Set Evolution |
    %|    One loop, alpha = 0     |
    %+----------------------------+
    fprintf('\tRemove Fragments by Level Set - Iteration # (refinement loop): %d\n\n', iter_outer + 1);
    alfa=0;
    phi = drlse_Denoise(phi, g, lambda, mu, alfa, epsilon, timestep, iter_inner, potentialFunction);


    varargout{1,1} = phi;
end

