%Sample Code
clear 
close all
keps = 1e-20;

d = 8;
m = 10*(d+1);
xlb = -10;
xub = 10;
ns = 1e-4;
xee = lhsdesign(m,d+1)'*(xub-xlb)+xlb;
% xee = rand(d+1,m)*(xub-xlb)+xlb;
fun = @(sig)some_noisy_function2(10.^sig(1:d),10^sig(d+1));

yee = zeros(1,m);
for i=1:m
    yee(i) = fun(xee(:,i));
end

xe=xee;
ye=yee;

% k = @(x1,x2,s)s(1)^2*exp(-1/(2*s(2)^2)*norm_matrix(x1,x2));  % squared exp
% k = @(x1,x2,s)s(1)^2*(1+sqrt(3)/s(2)^2*sqrt(norm_matrix(x1,x2))).*exp(-sqrt(3)/s(2)^2*sqrt(norm_matrix(x1,x2)));  % % Matern Class nu=3/2
k = @(x1,x2,s)s(1)*(1+sqrt(5)/s(2)*sqrt(norm_matrix(x1,x2))+5/3*norm_matrix(x1,x2)/s(2)^2).*exp(-sqrt(5)/s(2)*sqrt(norm_matrix(x1,x2)));  % % Matern Class nu=5/2
negml = @(sig)0.5*ye*inv(k(xe,xe,sig(1:2))+(ns^2)*eye(size(xe,2)))*ye' + 0.5*log((det(k(xe,xe,sig(1:2))+ns^2*eye(size(xe,2)))));
bestsol = cmaes(negml, 2, 1, 20, 30);
sig=bestsol.Position;
sig=[sig,ns];

inin=zeros(2,1);
for i=1:1000
    ymin = min(ye);
    Ke = k(xe,xe,sig(1:2))+(sig(3)^2)*eye(size(xe,2));
    invKe = inv(Ke);
    yh = @(x)k(xe,x,sig(1:2))*invKe*ye';
    s = @(x)sqrt(max(k(x,x,sig(1:2))-k(xe,x,sig(1:2))*invKe*k(xe,x,sig(1:2))',keps));
    EI = @(x)-1*((ymin-yh(x')) * cdf('Normal',(ymin-yh(x'))/s(x'),0,1) + s(x') * pdf('Normal',(ymin-yh(x'))/s(x'),0,1));
    
    bestsol = cmaes(EI, 9, xlb, xub, 80);
    xmin = bestsol.Position';
    xe = [xe,xmin];
    ym = fun(xmin);    
    ye = [ye,ym];
    
    i   
    
end

[ymin, n]=min(ye);
xmin=xe(:,n)