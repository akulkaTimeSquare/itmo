temp_out = imread("log_edited.png");
img_out = double(temp_out)/255;
abs_out = exp(img_out*m) - 1;
fourier_out = abs_out.*exp(1i.*an);
img_prepared = ifft2(ifftshift(fourier_out));
imwrite(ul, "edited.png");
imshow(img_prepared, []);