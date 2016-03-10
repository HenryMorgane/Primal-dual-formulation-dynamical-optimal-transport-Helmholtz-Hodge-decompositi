function [a,b]=proxJ2D(alpha,beta)
%
% proxJ3D - proximal operator of J in 2 dimensions
% projection onto the paraboloid K={alpha+(beta.^2)/2<=0}
%
%
% J(alpha,beta) = sum_i ||beta^i||^2/(alpha^i) + \Chi_{alpha^i>epsilon}
%
%

 a=alpha*0; b = beta*0; nb = beta*0;
 p13 = @(x) sign(x).*(abs(x).^(1./3));
 
 sgn = alpha+(beta.^2)/2;
 id = find(sgn<=0); % we are in K
 if (isempty(id)==0)
     a(id)=alpha(id); b(id)=beta(id);
 end
 
 % Reduction to a depressed cubic
 p = 2*(alpha+1); 
 q = -2*abs(beta);
 DELTA = q.^2+4./27*p.^3;
 ip = find((sgn>0)&(DELTA>0));
 if (isempty(ip)==0)

 	nb(ip) = p13((-q(ip)+sqrt(DELTA(ip)))/2)+p13((-q(ip)-sqrt(DELTA(ip)))/2);
 	a(ip) = -(nb(ip).^2)/2;
    b(ip) = sign(beta(ip)).*nb(ip);
 end
 
 im = find((sgn>0)&(DELTA<=0));
 if (isempty(im)==0)
 % we are outside of the paraboloid so the projection has a unique solution
 % what we calculate can have three solutions because it's the intersection
 % with a perpendicular (it can happen in three places)
 % two of the tree are on the negative part of the paraboloid so we consider
 % the positive solution
 
%      nb(im) = max(2*sqrt(-p(im)/3).*cos(1./3*acos(-(q(im)/2).*sqrt(-27./(p(im).^3)))),...
%               max(2*sqrt(-p(im)/3).*cos(1./3*acos(-(q(im)/2).*sqrt(-27./(p(im).^3)))+2*pi/3),...
%               2*sqrt(-p(im)/3).*cos(1./3*acos(-q(im)/2.*sqrt(-27./(p(im).^3)))+4*pi/3)));
% Faster below
     nb(im) = 2*sqrt(-p(im)/3).*cos(1./3*acos(-(q(im)/2).*sqrt(-27./(p(im).^3))));
     imm = find((sgn>0)&(nb<0)&(DELTA<=0));
     if (isempty(imm)==0)
         nb(imm)=2*sqrt(-p(imm)/3).*cos(1./3*acos(-(q(imm)/2).*sqrt(-27./(p(imm).^3)))+2*pi/3);
         immm = find((sgn>0)&(nb<0)&(DELTA<=0));
         if (isempty(immm)==0)
             nb(immm) = 2*sqrt(-p(immm)/3).*cos(1./3*acos(-q(immm)/2.*sqrt(-27./(p(immm).^3)))+4*pi/3);
         end
     end
     a(im) = -(nb(im).^2)/2;
  	 b(im) = sign(beta(im)).*nb(im);
     id = find((nb(im)<0),1);
     if (isempty(id)==0) 
        fprintf('Anormal case in projK2d.m\n');
     end
 end

 
