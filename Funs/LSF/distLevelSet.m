function [ varargout ] = distLevelSet( im, phi1, phi_mask_another, h_c, ...
                                       iter_inner, iter_outer, kappa, chi)
%   Detailed explanation goes here

    Img=double(im(:,:,1));
    %==================================
    %       Parameter setting
    %==================================
    timestep=5;  % time step
    mu=0.2/timestep;  % coefficient of the distance regularization term R(phi)
    epsilon=1.5; % papramater that specifies the width of the DiracDelta function
 
    %==================================
    %        Edge indicator
    %==================================
    sigma=1.5;     % scale parameter in Gaussian kernel
    G=fspecial('gaussian',15,sigma);
    Img_smooth=conv2(Img,G,'same');  % smooth image by Gaussiin convolution
    [Ix,Iy]=gradient(Img_smooth);
    f=Ix.^2+Iy.^2;
    g=1./(1+f);  % edge indicator function.
 
    %==================================
    % Initialize LSF as binary step 
    % function
    %==================================
    c0=2;
    initialLSF=c0*ones(size(Img));
    % generate the initial region R0 as a rectangle
    initialLSF(phi1 == 1)=-c0;  
    phi=initialLSF;
    
    
    %==================================
    %   Choose potential function
    %==================================
    potential=2;  
    if potential ==1
        potentialFunction = 'single-well';  % use single well potential p1(s)=0.5*(s-1)^2, which is good for region-based model 
    elseif potential == 2
        potentialFunction = 'double-well';  % use double-well potential in Eq. (16), which is good for both edge and region based models
    else
        potentialFunction = 'double-well';  % default choice of potential function
    end

    %==================================
    %   Level set evolution - outer
    %==================================
    for n=1:iter_outer
        
        fprintf('\t\tLSF Evolution - iteration #: %d\n', n);
        
        %==================================
        %            Update phi
        %==================================
        phi = dist2_drlse_edge(phi, phi_mask_another, g, h_c, ...
                               kappa, mu, epsilon, timestep, iter_inner, chi, potentialFunction); 
    end
    disp(' ');
    
    varargout{1,1} = phi;
end
