function LR=SVT(temp,ui,alpha)

   
    [U, sigma, V] = svd(temp, 'econ');
    diagS = diag(sigma);
    svp = length(find(diagS > alpha / ui));
    if svp >= 1
        diagS = diagS(1 : svp) - alpha / ui;
    else
        svp = 1;
        diagS = 0;
    end
    LR = U(:, 1 : svp) * diag(diagS) * V(:, 1:svp)';