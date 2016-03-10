function [dxu, dyu, dzu] = Grad(u)

d = size(u);

dxu = diff(u,[],1)*d(1);
dyu = diff(u,[],2)*d(2);
dzu = diff(u,[],3)*d(3);
