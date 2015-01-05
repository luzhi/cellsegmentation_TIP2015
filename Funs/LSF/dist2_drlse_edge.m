function phi_i = dist2_drlse_edge(phi_0, phi_j, g, h_c, ...
                                kappa, mu, epsilon, timestep, iter, chi, potentialFunction)

%==========================================
% Assign 'phi_0' to phi in inner iteration
%==========================================
phi_i=phi_0;
phi_j = ~phi_j;

% % %==================================
% % %  Gradient of edge indicator 'g'
% % %==================================
% % [vx, vy]=gradient(g);

for k=1:iter
    phi_i=NeumannBoundCond(phi_i);
    [phi_i_x,phi_i_y]=gradient(phi_i);
    s=sqrt(phi_i_x.^2 + phi_i_y.^2);
    smallNumber=1e-10;  
    Nx=phi_i_x./(s+smallNumber); % add a small positive number to avoid division by zero
    Ny=phi_i_y./(s+smallNumber);
    curvature=div(Nx,Ny);
    if strcmp(potentialFunction,'single-well')
        distRegTerm = 4*del2(phi_i)-curvature;  % compute distance regularization term in equation (13) with the single-well potential p1.
    elseif strcmp(potentialFunction,'double-well');
        distRegTerm=distReg_p2(phi_i);  % compute the distance regularization term in eqaution (13) with the double-well potential p2.
    else
        disp('Error: Wrong choice of potential function. Please input the string "single-well" or "double-well" in the drlse_edge function.');
    end        
    
    
    diracPhi_i=Dirac(phi_i,epsilon);

% %     %==================================
% %     %           Edge Term
% %     %==================================
% %     edgeTerm = diracPhi_i.*(vx.*Nx+vy.*Ny) + diracPhi_i.*g.*curvature;
    
    %==================================
    %   Dist map based length Term
    %==================================
%     w_dist_pair_g = g .* h_c .* clumpMask;
    w_dist_pair_g = g .* h_c;
    [vxd, vyd] = gradient(w_dist_pair_g);
%     distTerm = diracPhi_i.*(vxd.*Nx+vyd.*Ny) + diracPhi_i .* h_c .* clumpMask .* curvature;
    distTerm = diracPhi_i.*(vxd.*Nx+vyd.*Ny) + diracPhi_i .* h_c .* g .* curvature;        
   
    %+-------------------------+
    %|  Overlapping Area Term  |
    %+-------------------------+
    if ~isempty(intersect(find(~im2bw(phi_i) == 1), find(~phi_j == 1)))
        Heaviside_phi_j = Heaviside(phi_j);
        OverAreaTerm = h_c .* g .* diracPhi_i .* Heaviside_phi_j;
    else
        OverAreaTerm = 0;
    end
    
    %==================================
    %           Update phi
    %==================================
% %     phi_i = phi_i + timestep*(mu * distRegTerm + ...      % Regularization term
% %                               lambda * edgeTerm + ...   % Edge term                      
% %                               kappa * distTerm + ...    % Dist map length term
% %                               cha * OverAreaTerm);      % Overlapping area term
     phi_i = phi_i + timestep*(mu * distRegTerm + ...      % Regularization term                     
                               kappa * distTerm + ...    % Dist map length term
                               chi * OverAreaTerm);      % Overlapping area term
end

%==================================
%       Internal functions
%==================================

%+--------------------------+
%|    Heaviside function    |
%+--------------------------+
function f = Heaviside(phi)
% covert distance function phi to a binary Heaviside mask
% f = ~im2bw(phi);
phi_logical_tmp = false(size(phi));
phi_logical_tmp(phi <= 0) = 1;
f = phi_logical_tmp;

%==================================
%           Area of phi
%==================================
function f = phiArea(phi)
% count the number of pixels inside phi
phiBw = ~im2bw(phi);
f = length(find(phiBw == 1));



function f = distReg_p2(phi)
% compute the distance regularization term with the double-well potential p2 in eqaution (16)
[phi_x,phi_y]=gradient(phi);
s=sqrt(phi_x.^2 + phi_y.^2);
a=(s>=0) & (s<=1);
b=(s>1);
ps=a.*sin(2*pi*s)/(2*pi)+b.*(s-1);  % compute first order derivative of the double-well potential p2 in eqaution (16)
dps=((ps~=0).*ps+(ps==0))./((s~=0).*s+(s==0));  % compute d_p(s)=p'(s)/s in equation (10). As s-->0, we have d_p(s)-->1 according to equation (18)
f = div(dps.*phi_x - phi_x, dps.*phi_y - phi_y) + 4*del2(phi);  

function f = div(nx,ny)
[nxx,junk]=gradient(nx);  
[junk,nyy]=gradient(ny);
f=nxx+nyy;

function f = Dirac(x, sigma)
f=(1/2/sigma)*(1+cos(pi*x/sigma));
b = (x<=sigma) & (x>=-sigma);
f = f.*b;

function g = NeumannBoundCond(f)
% Make a function satisfy Neumann boundary condition
[nrow,ncol] = size(f);
g = f;
g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);  