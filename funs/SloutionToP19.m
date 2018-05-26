%
%  min  1/2 sum_v|| s - qv||^2
%  s.t. s>=0, 1's=1
function [x, ft] = SloutionToP19(q0, m)

if nargin < 2
    m = 1;
end;
ft=1;
n = length(q0);
p0 = sum(q0,1)/m-mean(sum(q0,1))/m + 1/n;
vmin = min(p0);
if vmin < 0
    f = 1;
    lambda_m = 0;
    while abs(f) > 10^-10
        v1 = lambda_m-p0;
        posidx = v1>0;
        npos = sum(posidx);
        g = npos/n-1;
        if 0 == g
            g = eps;
        end;
        f = sum(v1(posidx))/n - lambda_m;
        lambda_m = lambda_m - f/g;
        ft=ft+1;
        if ft > 100
            x = max(-v1,0);
            break;
        end;
    end;
    x = max(-v1,0);
else
    x = p0;
end;
