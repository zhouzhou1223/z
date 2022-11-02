% q: irf->srf
% q1,q2,q3: irf->ssrf
%
function q = Qcomination(fval1, q1, fval2, q2, fval3, q3)

R1_s = [0.999991953964000,-0.003855453067860,0.001107921250810;
    -0.002875276132160,-0.496285685373000,0.868154508875000;
    -0.002797283507320,-0.868150709252000,-0.496292777733000];    % ssrf,1->srf
R2_s = [0.999868439135000,0.015726793513000,-0.003971446564830;
    0.016149312081100,-0.942268716879000,0.334468032720000;
    0.001517939828470,-0.334488165946000,-0.942398728087000];
R3_s = [0.011846242780200,-0.769183928773000,0.638917639645000;
    -0.491411293086000,0.551999304112000,0.673655482637000;
    -0.870847063243000,-0.321951629871000,-0.371446551289000]; 

q1 = Qcontinues(fval1,q1);
q2 = Qcontinues(fval2,q2);
q3 = Qcontinues(fval3,q3);

lamda1 = 0.3; 
lamda2 = 0.3;
lamda3 = 0.3; % max(res, 0.3) the scaling factor for the larger uncertainty in rotation wrt. the boresight
k = 10; % kappa
M2 = R2_s*[lamda2^2,0,0;0,lamda2^2,0;0,0,lamda2^2*k^2]^(-1)*R2_s';
M3 = R3_s*[lamda3^2,0,0;0,lamda3^2,0;0,0,lamda3^2*k^2]^(-1)*R3_s';
M1 = R1_s*[lamda1^2,0,0;0,lamda1^2,0;0,0,lamda1^2*k^2]^(-1)*R1_s';

q1_i_s = quatmultiply(q1, RtoQ(R1_s)); % irf->ssrf, ssrf-> srf == irf->srf
q2_i_s = quatmultiply(q2, RtoQ(R2_s));
q3_i_s = quatmultiply(q3, RtoQ(R3_s)); 

N = max([length(q1),length(q2),length(q3)]); 

q = zeros(N,4);  % q crf->irf

for n=1:N
    if fval1(n,:)==1 && fval2(n,:)==1 && fval3(n,:)==0
        q12 = quatmultiply(quatinv(q1_i_s(n,:)), q2_i_s(n,:));
        d12 = 2*sign(q12(:,1))*[q12(:,2),q12(:,3),q12(:,4)]';
        s = 1/2*(M1+M2)^(-1)*(M2*d12);
        q(n,:) = quatmultiply(q1_i_s(n,:), [sqrt(abs(1-norm(s)^2)), s']);
    end
    if fval1(n,:)==1 && fval2(n,:)==0 && fval3(n,:)==1
        q13 = quatmultiply(quatinv(q1_i_s(n,:)), q3_i_s(n,:));
        d13 = 2*sign(q13(:,1))*[q13(:,2),q13(:,3),q13(:,4)]';
        s = 1/2*(M1+M3)^(-1)*(M3*d13);
        q(n,:) = quatmultiply(q1_i_s(n,:), [sqrt(abs(1-norm(s)^2)), s']);
    end
    if fval1(n,:)==0 && fval2(n,:)==1 && fval3(n,:)==1
        q23 = quatmultiply(quatinv(q2_i_s(n,:)), q3_i_s(n,:));
        d23 = 2*sign(q23(:,1))*[q23(:,2),q23(:,3),q23(:,4)]';
        s = 1/2*(M2+M3)^(-1)*(M3*d23);
        q(n,:) = quatmultiply(q2_i_s(n,:), [sqrt(abs(1-norm(s)^2)), s']);
    end
    if fval1(n,:)==1 && fval2(n,:)==1 && fval3(n,:)==1
        q12 = quatmultiply(quatinv(q1_i_s(n,:)), q2_i_s(n,:));
        d12 = 2*sign(q12(:,1))*[q12(:,2),q12(:,3),q12(:,4)]';
        q13 = quatmultiply(quatinv(q1_i_s(n,:)), q3_i_s(n,:));
        d13 = 2*sign(q13(:,1))*[q13(:,2),q13(:,3),q13(:,4)]';
        s = 1/2*(M1+M2+M3)^(-1)*(M2*d12+M3*d13);
        q(n,:) = quatmultiply(q1_i_s(n,:), [sqrt(abs(1-norm(s)^2)), s']);
    end
end
q = quatnormalize(q);
end