function dx = derivative(x)
%DERIVATIVE First-difference derivative with same-length output.
%
%   dx = derivative(x)
%
%   Input:
%     x  - signal vector
%   Output:
%     dx - approximate derivative (same length as x)

    arguments
        x (:,1) double
    end
    
    dx = diff(x);
    dx(end+1,1) = dx(end);   % pad to keep same length
end
