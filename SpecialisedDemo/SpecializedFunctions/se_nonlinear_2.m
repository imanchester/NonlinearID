function [se,ys,xs] = se_nonlinear_2(u,y,x1,e,f,g,ec,fc,gc,varargin)
% Compute simulation error for generic nonlinear model
% 30/11/2016 ~ include inputs in g.

options = optimset('Display','off');

[nx,~] = size(x1);
[ny,T] = size(y);

xs = zeros(nx,T);
ys = zeros(ny,T);

xt = x1;

for t = 1:T
 
% Store state, output:  
    xs(:,t) = xt;  
    ys(:,t) = g(xt,u(:,t),gc);
    
% Use a solver for nonlinear dynamics:
     [xt,~] = fsolve(@(x)(e(x,ec) -  f(xt,u(:,t),fc)),xt,options);
     
end 

% Update: if we supply an additional argument, treat it as a scaling factor
if length(varargin) == 1
    osf = varargin{1};
    ys = osf*ys;
end

tmp = y-ys;
tmp = tmp(:);
se = tmp'*tmp/(y(:)'*y(:));


end

