function [sp] = sigma_points(x, P, alpha)

%Generates sigma points given a mean state and a covariance matrix.

n = length(x);
k = 0;
lambda = alpha^2*(n+k)-n;
sq_mat = chol((n+lambda)*P, 'lower');
sp = zeros(n, 2*n);

for i = 1:(2*n)
    if i <= n
        if size(x) ~= size(sq_mat(:, i))
            a = 2;
        end
        sp(:, i) = x + sq_mat(:, i);
    else
        sp(:, i) = x - sq_mat(:, i-n);
    end
end
sp = [x, sp];
end
