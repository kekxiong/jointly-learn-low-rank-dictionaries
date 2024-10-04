function D2 = step2(H1, A1, gamma, u ,rho)

maxlamda=1e+6;

stop = false;

max_iter = 1000;

thre = 1e-6;


[atoms ,~] = size(A1);




D2 = H1(:,1:atoms);

y4 = zeros(size(D2));
O = zeros(size(D2));



for j=1:max_iter

    

    

    D2 = (H1 * A1'+ u * O + y4)/(A1 * A1'+ u * eye(size(A1 * A1')));
    

    O = SVT((D2 - y4/u),u ,gamma);
 

    

    y4 = y4 + min(u,maxlamda) * (O - D2);

    u=u * rho;
    
    %º∆À„rmse
    



    res = norm(O - D2, 'fro');
    
   

    if res < thre 
        stop = true;
        break;
    end
end