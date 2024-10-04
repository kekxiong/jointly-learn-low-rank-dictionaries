function DD = Deg_Dist(ref,tar)

%Degree of Distortion 
%���������ͼ����ԭʼͼ�����������𣬷�Ӧ�ع�ͼ��ʧ��̶ȸߵ�
% INPUT
%   ref : reference HS data (rows,cols,bands)
%   tar : target HS data (rows,cols,bands)
%
%   
%  [1]Zhu X X, Bamler R. A sparse image fusion algorithm with application to pan-sharpening[J].
%            IEEE transactions on geoscience and remote sensing, 2012, 51(5): 2827-2836.

%__________________________________________________________________________

[rows,cols,bands] = size(ref);
out = zeros(1,bands);
 for i = 1:bands
        tar_tmp = tar(:,:,i);
        ref_tmp = ref(:,:,i);
        dd_per = sum(sum(abs(ref_tmp - tar_tmp)))/(rows * cols);
        out(1,i) = dd_per;
 end
 
DD = mean(out);   