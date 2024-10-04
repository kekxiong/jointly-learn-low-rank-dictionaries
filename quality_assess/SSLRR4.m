function [D_h ,D_m] = SSLRR4( H ,M ,alpha1,beta1,alpha2,beta2,atoms)

% Xue J, Zhao Y Q, Bu Y, et al. 
% Spatial-spectral structured sparse low-rank representation for hyperspectral image super-resolution[J]. 
% IEEE Transactions on Image Processing, 2021, 30: 3084-3097.
% X = D_m *A_m +E_m;
%  min ||A||_* +alpha*||A||_1 +beta1*||E_h||_2,1 +beta2*||E_m||_2,1
% s.t. H = D_h *A +E_h ,M = D_m *A +
%___________________________________________________________________________




%参数设置
u_max = 1e+6 ;
u = 0.001 ;
lamda = 1.5;
thre = 1e-6;
max_iter = 100;
stop = false;
iter = 1;
delta=1e-6;

D_m=M(:,1:atoms);
D_h=H(:,1:atoms);
[bandms, pixels] = size(M);
[bandhs, ~] = size(H);
[~ , atoms] = size(D_m);







%变量初始化
A_h = zeros(atoms , pixels); %稀疏系数
A_m = zeros(atoms , pixels);
J = A_h; K = A_m;                     %辅助变量
P1 = zeros(size(A_h)); P2 = P1; P3 = P1; P4 = P1;P5 = P1;
[U1 ,S1,V1] = svd(A_h ,'econ');[U2 ,S2,V2] = svd(A_m ,'econ');
Uh = U1 ; Vh = V1;
Um = U2 ; Vm = V2;
while ~stop && iter < max_iter+1
    
    
    %更新U_h 和U_m
    iter
    temph1 = (A_h + P2/u) * Vh;
    Uh = sign(temph1).* max(abs(temph1)- beta1/u ,0);
    tempm1 = (A_m + P4/u) * Vm;
    Um = sign(tempm1).* max(abs(tempm1)- beta2/u ,0);
    
    %更新V
    temph2 = (A_h + P2/u)' * Uh;
    [Bh ,Dh ,Ch] = svd(temph2 ,'econ');
    Vh = Bh * Ch';
    tempm2 = (A_m + P4/u)' * Um;
    [Bm ,Dm ,Cm] = svd(tempm2 ,'econ');
    Vm = Bm * Cm';    
    
    %更新J
    temph3 = A_h +P1/u; J = SVT(temph3,u,alpha1);
    tempm3 = A_m +P3/u; K = SVT(tempm3,u,alpha2);
    
    %更新A_m

    h1 = D_h' * D_h +3* u* eye(atoms);
    h2 = D_h' * H+ u* (J +Uh *Vh'+ A_m)- P1 -P2 -P5;
    h3 = pinv(h1) * ones(atoms,1)*pinv(ones(1, atoms)* pinv(h1)* ones(atoms,1));
    A_h = pinv(h1) * (h2) -h3 * (ones(1,atoms)*pinv(h1)* h2 - 1);
 
    m1 = D_m' * D_m +3* u* eye(atoms);
    m2 = D_m' * M+ u* (K +Um *Vm'+ A_h)- P3 -P4 +P5 ;
    m3 = pinv(m1) * ones(atoms,1)*pinv(ones(1, atoms)* pinv(m1)* ones(atoms,1));
    A_m = pinv(m1) * (m2) -m3 * (ones(1,atoms)*pinv(m1)* m2 - 1);    
    
    %更新D_h 和D_m
%     D_h = (H*A_h') *pinv(A_m * A_h');
%     D_m = (M*A_m') *pinv(A_m * A_m');
    for j=1:atoms
        phi_h(j)=A_h(j,:)*A_h(j,:)'+delta;
        phi_m(j)=A_m(j,:)*A_m(j,:)'+delta;   
    end
    D_h=D_h+(H*A_h')/(phi_h);
    D_m=D_m+(M*A_m')/(phi_m);
    
    
    rmse_h=(sum(sum((H-D_h*A_h).^2))/(220*4350)).^0.5
    rmse_m=(sum(sum((M-D_m*A_m).^2))/(10*4350)).^0.5
    
    %更新拉格朗日乘数矩阵P1，P2，P3,P4

    P1= P1 + u* (A_h - J);
    P2= P2 + u* (A_h - Uh * Vh');
    P3= P3 + u* (A_m - K);
    P4= P4 + u* (A_m - Um * Vm');
    P5= P5 + u* (A_h - A_m);
    
    u = min(lamda * u ,u_max);
    
    iter=iter+1 ;   
    res1=norm(A_h - J, 'fro');
    res2=norm(A_h - Uh * Vh', 'fro');
    res3=norm(A_m - K, 'fro');
    res4=norm(A_m - Um * Vm', 'fro');
    res5=norm(A_h - A_m, 'fro');
    
    res = max([res1,res2, res3, res4 ,res5])
    if res < thre 
        stop = true;
        break;
    end
end

end