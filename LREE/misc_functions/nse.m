function se = nse(yd,ys)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

tmp = yd-ys;
tmp = tmp(:);
se = tmp'*tmp/(yd*yd');


end

