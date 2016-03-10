function a = proxJ3D(alpha)
%
% proxJ3D - proximal operator of J in 3 dimensions
%
% size of alpha = (rho,m) : (N,P,Q,3) 
% size of a : (N,P,Q,3) 
%
% J(rho,m) = sum_i ||m^i||^2/(rho^i) + \Chi_{rho^i>epsilon}
%

[a(:,:,:,3),~]=proxJ2D(alpha(:,:,:,3),sqrt(alpha(:,:,:,1).^2+alpha(:,:,:,2).^2));
 mu = alpha(:,:,:,3)-a(:,:,:,3);
 a(:,:,:,1) = alpha(:,:,:,1)./(1+mu);
 a(:,:,:,2) = alpha(:,:,:,2)./(1+mu);
 