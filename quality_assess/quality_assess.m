function [RMSE,PNSR,SAM,SSIM,ERGAS,UIQI,CC,DD]=quality_assess(x,y)

% y为参考图像
% x为重构的图像

%rows,columns,bands
[rows,columns,bands]=size(x);



% RMSE 
aux = sum(sum((x - y).^2, 1), 2)/(rows*columns);
rmse_per_band = sqrt(aux);
RMSE = sqrt(sum(aux, 3)/bands);

%PNSR
%PNSR = PSNR(x,y);
PNSR = csnr(y,x,0,0);



% SAM
sam = SpectAngMapper(y, x );
SAM = sam*180/pi;

%SSIM
M = zeros(1, bands);
for i = 1 : bands
    M(1, i) = ssim_index(255*mat2gray(x(:, :, i)), 255 * mat2gray(y(:, :, i)));
end
SSIM = mean(M);

%ERGAS
ERGAS = ErrRelGlobAdimSyn(x,y);

%UIQI
q_band = zeros(1, bands);
for idx1=1:bands
    q_band(idx1)=img_qi(y(:,:,idx1), x(:,:,idx1), 32);
end
UIQI = mean(q_band);

%Correlation Coefficient
CC = CoCo(y,x);

%DD(Degree of Distortion)
DD = Deg_Dist(y, x);


