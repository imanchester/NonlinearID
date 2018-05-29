function res = normalizeData(data,varargin)

% Do you want the data inside [-1,1] or to have unit variance?

if length(varargin) == 1
    ops = varargin{1};
else
    ops = [];
end

if isfield(ops,'u')
    norm_input = ops.u;
else
    norm_input = 1;
end

if isfield(ops,'y')
    norm_ouyput = ops.y;
else
    norm_ouyput = 1;
end

% Normalize inputs
u = data.u;
u0 = mean(u);
u_ = u - u0;
usf = diag(max(abs(u_),[],2)); 

u = usf\u_;

% Normalize outputs
y = data.y;
y0 = mean(y);
y_ = y - y0;
ysf = diag(max(abs(y_),[],2)); 

y = ysf\y_;

% Structure for correcting the data:
normvals.u0 = u0;
normvals.y0 = y0;
normvals.usf = usf;
normvals.ysf = ysf;

% Finally, return 
data.u = u;
data.y = y;

res.data = data;
res.normvals = normvals;


















