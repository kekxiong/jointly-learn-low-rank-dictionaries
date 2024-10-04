function [D1 ,A1 ,A2] = first_step(M1, M2, atoms, alpha,beta,lamda1,lamda2,u,rho)

maxlamda=1e+6;

stop = false;

max_iter = 1000;

thre = 1e-6;


[bandm,n1] = size(M1);
[~,n2] = size(M2);




D1 = M1(:,1:atoms);


y1 = zeros(atoms,n1); 
y2 = zeros(atoms,n2);
y3 = zeros(bandm, atoms);


J = zeros(size(y1));
K = zeros(size(y2));

L=zeros(size(D1));



for i=1:max_iter

    

    t1 = D1'* D1 + u * eye(size(D1'* D1));
    t2 = D1'* M1 + y1 + u * J;
    t3 = pinv(t1) * ones(atoms, 1) * pinv(ones(1,atoms) * pinv(t1) * ones(atoms, 1));
    A1 = pinv(t1) * t2 - t3 * (ones(1, atoms) * pinv(t1) * t2 - 1);

 
    s1 = alpha * D1'* D1 + u * eye(size(D1'* D1));
    s2 = alpha * D1'* M2 + y2 + u * K;
    s3 = pinv(s1) * ones(atoms, 1) * pinv(ones(1,atoms) * pinv(s1) * ones(atoms, 1));
    A2 = pinv(s1) * s2 - s3 * (ones(1, atoms) * pinv(s1) * s2 - 1); 
    

    D1 = (M1 * A1'+ alpha * M2 * A2' + u* L + y3)/(A1 * A1'+ alpha * A2 * A2'+ u*eye(size(A1 * A1')));
    

    L = SVT((D1 - y3/u),u ,beta);
  

    J = max(abs(A1 - y1/u)- lamda1/u, 0).*sign(A1 - y1/u);
    K = max(abs(A2 - y2/u)- lamda2/u, 0).*sign(A2 - y2/u);



    
    y1 = y1 + min(u,maxlamda) * (J - A1);
    y2 = y2 + min(u,maxlamda) * (K - A2);
    y3 = y3 + min(u,maxlamda) * (L - D1);

    u=u * rho;
    

    

    

    res1 = norm(J - A1, 'fro');
    res2=norm(K - A2, 'fro');
    res3=norm(L - D1, 'fro');
    

   

    if res1 < thre && res2 < thre && res3 < thre 
        stop = true;
        break;
    end
end


end