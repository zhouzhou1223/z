%
% making quaternions continuous
%
% INPUT： Q : Quaternion
%         fQ: flag for quaternions, 
%
% 
function Q = Qcontinues(fQ,Q)
N = size(Q, 1);
for n = 2 : N
    if fQ(n) == 0
        Q(n,:) = Q(n-1, :);
    else 
        if Q(n,:)*Q(n-1, :)' < 0
            Q(n,:) = -Q(n, :);
        end
    end
end
    