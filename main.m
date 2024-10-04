clear all ;
close all;
addpath('funs');
addpath('quality_assess');
addpath('Data');



load MSI_IndinePine.mat;
HSI = double(imread('19920612_AVIRIS_IndianPine_Site3.tif'));
MSI=data_MS_HR;

MSI = MSI / max(max(max(MSI)));
HSI = HSI / max(max(max(HSI)));

train_siz1 = 73;

M1 = MSI(:, 1 :train_siz1, :);
H1 = HSI(:,1 :train_siz1, :);
M2=MSI(:, train_siz1+1: end, :);
H2=HSI(:, train_siz1+1: end, :);

M1_2d = hyperConvert2d(M1);
H1_2d = hyperConvert2d(H1);
M2_2d = hyperConvert2d(M2);
H2_2d = hyperConvert2d(H2);

[h, w, bandhs] = size(H2);


%≤Œ ˝…Ë÷√
atoms = 100;
alpha = 3;
beta = 0.08;
lamda1 = 0.005;
lamda2 = 0.005;

u = 0.001;
rho = 1.4;
gamma = 0.05;







[D1 ,A1 ,A2] = first_step(M1_2d, M2_2d, atoms, alpha,beta,lamda1,lamda2,u,rho);

D2 = step2(H1_2d , A1, gamma, u ,rho);

rc = D2 * A2 ;

rc_HSI = hyperConvert3d(rc, h, w, bandhs);

[RMSE,PNSR,SAM,SSIM,ERGAS,UIQI,CC,DD]=quality_assess(rc_HSI,H2)




