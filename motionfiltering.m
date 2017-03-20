clear all
close all
clc
imtool close all

image_set = 1;

switch image_set
    case 1
        imstring = 'EnterExitCrossingPaths2cor%04d.jpg';
        im_start = 0;
        im_end = 484;
    case 2
        imstring = 'img01_%04d.jpg';
        im_start = 1;
        im_end = 1070;
    case 3
        imstring = 'advbgst1_21_%04d.jpg';
        im_start = 2;
        im_end = 354;
end

num_images = im_end-im_start+1;
for i = im_start:im_end
    images(:,:,i+1) = rgb2gray(imread(sprintf(imstring,i)));
end

filter = [-1 0 1];
filter_depth = (length(filter)-1)/2;
filtered_image = zeros(size(images(:,:,1),1),size(images(:,:,1),2),...
    num_images);
mask = filtered_image;

for i = 1+filter_depth:num_images-filter_depth
    filtered_values = zeros(size(images(:,:,1)));
    for l = 1:length(filter)
        filtered_values = filtered_values + ...
            double(images(:,:,i-filter_depth+l-1))*filter(l);
    end
    filtered_image(:,:,i) = abs(filtered_values);
end
threshold = 15;
mask(filtered_image >= threshold) = 1;
mask(filtered_image < threshold) = 0;

result = images;
result(mask==1) = 255;

implay(result)
%%
%simple 1D

filter = 0.5*[-1 0 1];
filter_depth = (length(filter)-1)/2;
filtered_image = zeros(size(images(:,:,1),1),size(images(:,:,1),2),...
    num_images);
mask = filtered_image;

for i = 1+filter_depth:num_images-filter_depth
    filtered_values = zeros(size(images(:,:,1)));
    for l = 1:length(filter)
        filtered_values = filtered_values + ...
            double(images(:,:,i-filter_depth+l-1))*filter(l);
    end
    filtered_image(:,:,i) = abs(filtered_values);
end
threshold = 15;
mask(filtered_image >= threshold) = 1;
mask(filtered_image < threshold) = 0;

result = images;
result(mask==1) = 255;

implay(result)

%%
% gaussian
tsigma = 1.2;
gfilter = fspecial('gaussian',[1,5],tsigma);
filter = gradient(gfilter);
filter_depth = (length(filter)-1)/2;
filtered_image = zeros(size(images(:,:,1),1),size(images(:,:,1),2),...
    num_images);
mask = filtered_image;

for i = 1+filter_depth:num_images-filter_depth
    filtered_values = zeros(size(images(:,:,1)));
    for l = 1:length(filter)
        filtered_values = filtered_values + ...
            double(images(:,:,i-filter_depth+l-1))*filter(l);
    end
    filtered_image(:,:,i) = abs(filtered_values);
end
threshold = 5;
mask(filtered_image >= threshold) = 1;
mask(filtered_image < threshold) = 0;

result = images;
result(mask==1) = 255;

implay(result)

%%
% higher tsigma

tsigma = 1;
gfilter = fspecial('gaussian',[1,round(5*tsigma)],tsigma);
filter = gradient(gfilter);
filter_depth = (length(filter)-1)/2;
filtered_image = zeros(size(images(:,:,1),1),size(images(:,:,1),2),...
    num_images);
mask = filtered_image;

for i = 1+filter_depth:num_images-filter_depth
    filtered_values = zeros(size(images(:,:,1)));
    for l = 1:length(filter)
        filtered_values = filtered_values + ...
            double(images(:,:,i-filter_depth+l-1))*filter(l);
    end
    filtered_image(:,:,i) = abs(filtered_values);
end
threshold = 5;
mask(filtered_image >= threshold) = 1;
mask(filtered_image < threshold) = 0;

result = images;
result(mask==1) = 255;

implay(result)

%%
%3x3 box filter
box_size = 3;
for i = 1:num_images
    images_box(:,:,i) = imfilter(images(:,:,i),...
        fspecial('average',[box_size,box_size]));
end

filter = 0.5*[-1 0 1];
filter_depth = (length(filter)-1)/2;
filtered_image = zeros(size(images_box(:,:,1),1),...
    size(images_box(:,:,1),2),num_images);
mask = filtered_image;

for i = 1+filter_depth:num_images-filter_depth
    filtered_values = zeros(size(images_box(:,:,1)));
    for l = 1:length(filter)
        filtered_values = filtered_values + ...
            double(images_box(:,:,i-filter_depth+l-1))*filter(l);
    end
    filtered_image(:,:,i) = abs(filtered_values);
end
threshold = 15;
mask(filtered_image >= threshold) = 1;
mask(filtered_image < threshold) = 0;

result = images_box;
result(mask==1) = 255;

implay(result)

%%
%5x5 box filter
box_size = 5;
for i = 1:num_images
    images_box(:,:,i) = imfilter(images(:,:,i),...
        fspecial('average',[box_size,box_size]));
end

filter = 0.5*[-1 0 1];
filter_depth = (length(filter)-1)/2;
filtered_image = zeros(size(images_box(:,:,1),1),...
    size(images_box(:,:,1),2),num_images);
mask = filtered_image;

for i = 1+filter_depth:num_images-filter_depth
    filtered_values = zeros(size(images_box(:,:,1)));
    for l = 1:length(filter)
        filtered_values = filtered_values + ...
            double(images_box(:,:,i-filter_depth+l-1))*filter(l);
    end
    filtered_image(:,:,i) = abs(filtered_values);
end
threshold = 15;
mask(filtered_image >= threshold) = 1;
mask(filtered_image < threshold) = 0;

result = images_box;
result(mask==1) = 255;

implay(result)

%%
%2D Gaussian
tsigma = 1;
for i = 1:num_images
    images_gaussian(:,:,i) = imfilter(images(:,:,i),...
        fspecial('gaussian',[5,5],tsigma));
end

filter = 0.5*[-1 0 1];
filter_depth = (length(filter)-1)/2;
filtered_image = zeros(size(images_gaussian(:,:,1),1),...
    size(images_gaussian(:,:,1),2),num_images);
mask = filtered_image;

for i = 1+filter_depth:num_images-filter_depth
    filtered_values = zeros(size(images_gaussian(:,:,1)));
    for l = 1:length(filter)
        filtered_values = filtered_values + ...
            double(images_gaussian(:,:,i-filter_depth+l-1))*filter(l);
    end
    filtered_image(:,:,i) = abs(filtered_values);
end
threshold = 15;
mask(filtered_image >= threshold) = 1;
mask(filtered_image < threshold) = 0;

result = images_gaussian;
result(mask==1) = 255;

implay(result)

%%
%threshold optimization

filter = 0.5*[-1 0 1];
filter_depth = (length(filter)-1)/2;
filtered_image = zeros(size(images(:,:,1),1),...
    size(images(:,:,1),2),num_images);
mask = filtered_image;

for i = 1+filter_depth:num_images-filter_depth
    filtered_values = zeros(size(images(:,:,1)));
    for l = 1:length(filter)
        filtered_values = filtered_values + ...
            double(images(:,:,i-filter_depth+l-1))*filter(l);
    end
    filtered_image(:,:,i) = abs(filtered_values);
end

sub_images = double(images(size(images,1)-20:size(images,1),1:20,:));
E = 1/num_images * sum(sub_images,3);
% E mean squared difference
E_msd = 0;
for i = 1:num_images
    E_msd = E_msd + (E - sub_images(:,:,i)).^2;
end
sigma = sqrt(1/(num_images - 1)*E_msd);
est_sigma = mean(mean(sigma));
threshold = 3*est_sigma;

mask(filtered_image >= threshold) = 1;
mask(filtered_image < threshold) = 0;

result = images;
result(mask==1) = 255;

implay(result)
