function image_out = normalized_image3(image)
max_value = max(max(image));
min_value = min(min(image));
rangle = repmat(max_value - min_value,[size(image,1),size(image,2),1]);
image_out = (image - repmat(min_value,[size(image,1),size(image,2),1]))./rangle;