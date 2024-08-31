% 此仿真对象参考 

close all;
clear;
clc;

%system matrices
A = [0.9499 -0.0017 0.0070 0.0005 0 0; 7.0373 0.4782 0.2243 0.1825 0.0079 0.0008;
    0.0094 0.0004 0.9927 -0.0148 0.0033 0.0002; 0.1035 0.0175 0.0031 0.8348 0.0704 0.0039;
    -0.1396 0.0073 0.0005 0.0122 0.9907 -0.0003;-0.6253 -0.00524 -0.0093 -3.4999 1.8816 -0.0169];
B = [-0.0003 0.0005; -0.0737 0.1149; 0.0003 0; 0.0106 0.0016; -0.0003 -0.0011; -0.0398 -0.0075];

alpha = 0.9; %



% Q = 100*eye(6);
Q = diag([1,1,1,0.1,0.1,0.1]);

R = 0.01*eye(2);
% digits(13);
% Q1 = [C'*Q*C -C'*Q*T; -T'*Q*C T'*Q*T];

N = 501;
Qn = 1*eye(6); % terminal state
Step = 300;


for i = 1:Step 
    u = randn(2,N-1)+sin(i);

    x(1:6,i,1) = randn(6,1);

    for k = 1:N-1
       x(1:6,i,k+1) = A * x(1:6,i,k) + B * u(:,k); % 状态量
       x(7:8,i,k) = u(:,k); % 控制量
    end
end

for k = 1: N-1
    for i = 1: Step 
        fai(k,i,:) = [x(1,i,k)^2  2*x(1,i,k)*x(2,i,k) 2*x(1,i,k)*x(3,i,k) 2*x(1,i,k)*x(4,i,k) 2*x(1,i,k)*x(5,i,k) 2*x(1,i,k)*x(6,i,k) 2*x(1,i,k)*x(7,i,k) 2*x(1,i,k)*x(8,i,k)... 
            x(2,i,k)^2 2*x(2,i,k)*x(3,i,k) 2*x(2,i,k)*x(4,i,k) 2*x(2,i,k)*x(5,i,k) 2*x(2,i,k)*x(6,i,k) 2*x(2,i,k)*x(7,i,k) 2*x(2,i,k)*x(8,i,k)...
            x(3,i,k)^2 2*x(3,i,k)*x(4,i,k) 2*x(3,i,k)*x(5,i,k) 2*x(3,i,k)*x(6,i,k) 2*x(3,i,k)*x(7,i,k) 2*x(3,i,k)*x(8,i,k)...
            x(4,i,k)^2 2*x(4,i,k)*x(5,i,k) 2*x(4,i,k)*x(6,i,k) 2*x(4,i,k)*x(7,i,k) 2*x(4,i,k)*x(8,i,k)...
            x(5,i,k)^2 2*x(5,i,k)*x(6,i,k) 2*x(5,i,k)*x(7,i,k) 2*x(5,i,k)*x(8,i,k)...
            x(6,i,k)^2 2*x(6,i,k)*x(7,i,k) 2*x(6,i,k)*x(8,i,k)...
            x(7,i,k)^2 2*x(7,i,k)*x(8,i,k)...
            x(8,i,k)^2];
    end
end

for k = N-1:-1:1
    
    if k == N-1
        for i = 1:Step 
            gamma(k,i) = (1/(alpha*alpha))*x(1:6,i,k+1)' * Q * x(1:6,i,k+1) + x(1:6,i,k)' * Q * x(1:6,i,k) + x(7:8,i,k)' * R * x(7:8,i,k);  
        end
    else
        for i = 1:Step
            gamma(k,i) =(1/(alpha*alpha))*x(1:6,i,k+1)' * P_bar(:,:,k+1) * x(1:6,i,k+1) + x(1:6,i,k)' * Q * x(1:6,i,k) + x(7:8,i,k)' * R * x(7:8,i,k);
        end
    end
    rank(squeeze(fai(k,:,:))'*squeeze(fai(k,:,:)))

    lamda(:,k) = inv(squeeze(fai(k,:,:))'*squeeze(fai(k,:,:))) * squeeze(fai(k,:,:))' * gamma(k,:)';
    H_xx(:,:,k) = [lamda(1,k) lamda(2,k) lamda(3,k) lamda(4,k) lamda(5,k) lamda(6,k);...
                   lamda(2,k) lamda(9,k) lamda(10,k) lamda(11,k) lamda(12,k) lamda(13,k);...
                   lamda(3,k) lamda(10,k) lamda(16,k) lamda(17,k) lamda(18,k) lamda(19,k);...
                   lamda(4,k) lamda(11,k) lamda(17,k) lamda(22,k) lamda(23,k) lamda(24,k);...
                   lamda(5,k) lamda(12,k) lamda(18,k) lamda(23,k) lamda(27,k) lamda(28,k);...
                   lamda(6,k) lamda(13,k) lamda(19,k) lamda(24,k) lamda(28,k) lamda(31,k)];

    H_xu(:,:,k) = [lamda(7,k) lamda(8,k);lamda(14,k) lamda(15,k);lamda(20,k) lamda(21,k);...
                   lamda(25,k) lamda(26,k);lamda(29,k) lamda(30,k);lamda(32,k) lamda(33,k)];
    H_ux(:,:,k) = H_xu(:,:,k)';

    H_uu(:,:,k) = [lamda(34,k) lamda(35,k);lamda(35,k) lamda(36,k)];

    K_bar(:,:,k) = -inv(H_uu(:,:,k))*H_ux(:,:,k);
    P_bar(:,:,k) = H_xx(:,:,k) + K_bar(:,:,k)'*H_ux(:,:,k)+H_xu(:,:,k)*K_bar(:,:,k)+K_bar(:,:,k)'*H_uu(:,:,k)*K_bar(:,:,k);
    NormFP_bar(k) = norm(P_bar(:,:,k),'fro');
end



for kk = N-1:-1:1
    if kk == N-1
        P(:,:,kk) = (alpha^2*Q+A'*Q*A-A'*Q*B*inv(alpha^2*R+B'*Q*B)*B'*Q*A)/alpha^2;
    else
        P(:,:,kk) = (alpha^2*Q+A'*P(:,:,kk+1)*A - A'*P(:,:,kk+1)*B*inv(alpha^2*R+B'*P(:,:,kk+1)*B)*B'*P(:,:,kk+1)*A)/alpha^2;
    end
    NormFP(kk) = norm(P(:,:,kk),'fro');
end

for kk = 1:N-1
    if kk == N-1
        K(:,:,kk) = inv(alpha^2*R + B'*Q*B)*B'*Q*A; %应该是K(N-1)
    else
        K(:,:,kk) = inv(alpha^2*R + B'*P(:,:,kk+1)*B)*B'*P(:,:,kk+1)*A; % 这里应该是P(k+1)，但是
    end
end
% x_bar(:,1) = [0.1; -0.1; -0.1; 0.1; 0.1; -0.1];
x_bar(:,1) = [1; -1; -1; -1; 1; 1];

for kk = 1:N-1
    u_bar(:,kk) = K_bar(:,:,k)*x_bar(:,kk);

	x_bar(:,kk+1)= A*x_bar(:,kk) + B*u_bar(:,kk);

end

t = 1:1:N;
figure(1)
plot(t,x_bar)

z = 1:1:N-1;
figure(2)
plot(z,u_bar)


% x1(:,1) = [0.1; -0.1; -0.1; 0.1; 0.1; -0.1];
x1(:,1) = [1; -1; -1; -1; 1; 1];

for kk = 1:N-1                                                             

    u1(:,kk) = -K(:,:,kk)*x1(:,kk);

    x1(:,kk+1)= A*x1(:,kk) + B*u1(:,kk);

end


figure(5)
plot(z,u1)
figure(6)
plot(t,x1)


