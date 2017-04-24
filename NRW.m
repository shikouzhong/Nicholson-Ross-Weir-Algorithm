function [er, ei, ur, ui] = NRW(newS11, newS21, f, sample, k, kc)
% 实现NRW算法
%   双向测量请两次调用

X = ((newS11).^2 - (newS21).^2 + 1) ./ (2 .* newS11);

R = zeros(length(f), 1);
for m = 1 : length(f)
    R(m) = X(m) + (X(m)^2 - 1).^(1/2);
    if abs(R(m)) >= 1
        R(m) = X(m) - (X(m)^2 - 1).^(1/2);
    end
end
T = (newS11 + newS21 - R) ./ (1 - (newS11 + newS21) .* R);

tempLog = log(1.0 ./ T);
detect = abs(diff(imag(tempLog)));
index = find(detect > pi);
for m = index + 1: length(tempLog)
    tempLog(m) = log(1.0 / T(m)) + j * 2 * pi * 1;;
end

Lambda_inv_2 = -(tempLog / (2 * pi * sample)).^2;
Lambda_inv = zeros(length(f), 1);
for m = 1 : length(f)
    Lambda_inv(m) = (Lambda_inv_2(m)).^(1/2);
    if real(Lambda_inv(m)) < 0
        Lambda_inv(m) = -(Lambda_inv_2(m)).^(1/2);
    end
end

mu =  2 * pi * Lambda_inv .* ((1 + R) ./ (1 - R)) .* ((k.^2 - kc.^2).^(-1/2));
epsilon = (2 * pi ./ k).^2 ./ mu .* (Lambda_inv_2 + (kc / (2 * pi)).^2);

er = real(epsilon);
ei = -imag(epsilon);
ur = real(mu);
ui = -imag(mu); 

end

