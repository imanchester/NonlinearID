function x_f = dhpf(x,a)
% Simple digital high pass filter.

x_f = zeros(size(x));
x_pf = 0;
x_p = 0;
for i = 1:size(x,2)
    x_f(:,i) = a*x_pf + a*(x(:,i) - x_p);
    x_pf = x_f(:,i);
    x_p = x(:,i);
end

end

