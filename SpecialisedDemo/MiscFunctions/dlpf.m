function x_f = dlpf(x,a)
% Simple digital low pass filter.

x_f = zeros(size(x));
x_p = 0;
for i = 1:size(x,2)
    x_f(:,i) = (1-a)*x_p + a*x(:,i);
    x_p = x_f(:,i);
end

end

