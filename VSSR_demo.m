%% VSSR 
%##########################################################################
%##########################  Copyright (c) 2020  ##########################
%#####  Created on 2020-07-24                               ###############
%#####  Author:Yanglei_Wu                                   ###############
%#####  Email:yanglei.wuu@gmail.com                         ###############
%#####  VSSR bias correction algorithm                      ###############
%#####  Version 1.0                                         ###############
%#####  Copyright (c) 2020 Yanglei.Wu. All rights reserved  ###############
%##########################################################################

% In this script, we also provide MICO bias correction algorithm to compare, %
% which proposed in the following paper:                                     %
%      C. Li, J.C. Gore, and C. Davatzikos,                                  %
%      "Multiplicative intrinsic component optimization (MICO) for MRI bias  %
%      field estimation and tissue segmentation", Magnetic Resonance Imaging % 
%      vol. 32 (7), pp. 913-923, 2014                                        %

clear all
close all

Dataset = 'T2_abd_FS'; 
                            
                            % 'T1_head'
                            % 'T2_head'
                            % 'T1_abd'
                            % 'T2_abd'
                            % 'T2_abd_bal'
                            % 'T2_abd_FS'

distribution      = 'Rice';  % noise distribution
                             %  'Gauss' --> Gaussian distribution
                             %  'Rice ' --> Rician Distribution
profile           = 'mp';    % BM4D parameter profile
                             %  'lc' --> low complexity
                             %  'np' --> normal profile
                             %  'mp' --> modified profile
                             % The modified profile is default in BM4D. For 
                             % details refer to the 2013 TIP paper.
do_wiener         = 1;       % Wiener filtering
                             %  1 --> enable Wiener filtering
                             %  0 --> disable Wiener filtering
verbose           = 1;       % verbose mode

estimate_sigma    = 0;       % enable sigma estimation
load(Dataset);
switch Dataset
    case 'T1_head'
        image = ifftshift(ifft2(ifftshift(kspace)),2);
        image = fliplr(imrotate(image,-90));
        body = image(:,:,4:5,:);
        image(:,:,4:5,:) = [];
        ft = fittype( 'poly44' );
    case 'T2_head'
        image = ifftshift(ifft2(ifftshift(kspace)));
        image = fliplr(imrotate(image,-90));
        body = image(:,:,9:10,:);
        image(:,:,9:10,:) = [];
        ft = fittype( 'a .* (5.*y.*y.*y -3.*y)./2 + b * (3.*y.*y -1).*x./2 + c .* (3.*y.*y -1)./2 + d .* y.*(3.*x.*x -1)./2 + e .* x.*y + f .* y + g .* (5.*x.*x.*x - 3.*x)./2 + h .* (3.*x.*x - 1)./2 + k .* x + l;', 'independent', {'x', 'y'}, 'dependent', 'z' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
    case 'T1_abd'
        image = ifftshift(ifft2(kspace),2);
        body = image(:,:,15:16,:);
        image(:,:,15:16,:) = [];
        ft = fittype( 'poly44' );
    case 'T2_abd'
        image = ifftshift(ifft2(kspace));
        body = image(:,:,7:8,:);
        image(:,:,7:8,:) = [];
        ft = fittype( 'poly44' );
    case 'T2_abd_bal'
        image = ifftshift(ifft2(kspace),2);
        body = image(:,:,15:16,:);
        image(:,:,15:16,:) = [];
        ft = fittype( 'poly44' );
    case 'T2_abd_FS'
        image = ifftshift(ifft2(kspace));
        body = image(:,:,7:8,:);
        image(:,:,7:8,:) = [];
        ft = fittype( 'poly44' );
end
surface = image;
im_comb_body = sqrt(sum(abs(body).^2,3));
im_comb_surface = sqrt(sum(abs(surface).^2,3));
sigma = std2(im_comb_body(10:20,10:20,1,1));
[im_comb_body, sigma_est] = bm4d(squeeze(im_comb_body), distribution, (~estimate_sigma)*sigma, profile, do_wiener, verbose);
im_comb_body = reshape(im_comb_body,[size(im_comb_body,1),size(im_comb_body,2),1,size(im_comb_body,3)]);

IIS = im_comb_body./im_comb_surface;
BW2 = zeros(size(im_comb_body));
se = strel('disk',5');
for i = 1:size(BW2,4)
    T2 = graythresh(im_comb_body(:,:,1,i));
    BW2(:,:,1,i) = im2bw(im_comb_body(:,:,1,i),T2);
    BW2(:,:,1,i) = imfill(BW2(:,:,1,i),'holes');
    BW2(:,:,1,i) = imclose(BW2(:,:,1,i),se);
    BW2(:,:,1,i) = imfill(BW2(:,:,1,i),'holes');
end
M = IIS .* BW2;
f44_res = zeros(size(M));
y = (1:size(M,1)) - round(size(M,1)/2);
x = (1:size(M,2)) - round(size(M,2)/2);
for i = 1:size(M,4)
    [xData, yData, zData] = prepareSurfaceData( x, y, M(:,:,1,i) );
    if exist('opts')
        f44 = fit( [xData, yData], zData, ft, opts);
    else
        f44 = fit( [xData, yData], zData, ft);
    end
    f44_res(:,:,1,i) = reshape(f44(xData,yData),[size(M,1),size(M,2)]);
end
IIS_C = im_comb_surface .* f44_res;
im_comb_body = normalized_image3(squeeze(im_comb_body));
im_comb_surface = normalized_image3(squeeze(im_comb_surface));
IIS_C = normalized_image3(squeeze(IIS_C));
figure('NumberTitle', 'off', 'Name', 'VSSR');
subplot(3,3,1),imshow(im_comb_body(:,:,1)),title('body coil image');
subplot(3,3,2),imshow(im_comb_surface(:,:,1)),title('surface coil image');
subplot(3,3,3),imshow(IIS_C(:,:,1)),title('corrected surface coil image');
subplot(3,3,4),imshow(im_comb_body(:,:,2));
subplot(3,3,5),imshow(im_comb_surface(:,:,2));
subplot(3,3,6),imshow(IIS_C(:,:,2));
subplot(3,3,7),imshow(im_comb_body(:,:,3));
subplot(3,3,8),imshow(im_comb_surface(:,:,3));
subplot(3,3,9),imshow(IIS_C(:,:,3));

%% MICO
% This code implements the MICO algorithm for joint segmentation and  bias field estimation 
% proposed in the following paper:      
%      C. Li, J.C. Gore, and C. Davatzikos, 
%      "Multiplicative intrinsic component optimization (MICO) for MRI bias field estimation and tissue segmentation", Magnetic Resonance
%      Imaging , vol. 32 (7), pp. 913-923, 2014
%
% All rights researved by Chunming Li
% E-mail: li_chunming@hotmail.com
% URL: http://imagecomputing.org/~cmli/
% Copyright (c) by Chunming Li
% Author: Chunming Li

clearvars -except im_comb_body im_comb_surface BW2

for i = 1:size(im_comb_surface,3)
    IIS_C(:,:,i) = Demo_MICO(im_comb_surface(:,:,i),BW2(:,:,i));
end
IIS_C = normalized_image3(squeeze(IIS_C));
figure('NumberTitle', 'off', 'Name', 'MICO');
subplot(3,3,1),imshow(im_comb_body(:,:,1)),title('body coil image');
subplot(3,3,2),imshow(im_comb_surface(:,:,1)),title('surface coil image');
subplot(3,3,3),imshow(IIS_C(:,:,1)),title('corrected surface coil image');
subplot(3,3,4),imshow(im_comb_body(:,:,2));
subplot(3,3,5),imshow(im_comb_surface(:,:,2));
subplot(3,3,6),imshow(IIS_C(:,:,2));
subplot(3,3,7),imshow(im_comb_body(:,:,3));
subplot(3,3,8),imshow(im_comb_surface(:,:,3));
subplot(3,3,9),imshow(IIS_C(:,:,3));