temp_initial = imread("toedit.png");
imshow(temp_initial, [])
img_initial = double(temp_initial)/255;
fourier = fftshift(fft2(img_initial));
ab = abs(fourier); % абсолютные значения
an = angle(fourier); % фазы
l = log(ab+1);
m = max(l(:));
ul = l/m; % нормировка логарифма
imwrite(ul, "log.png");
%imshow(ul, []);

temp_out = imread("log_edited.png");
img_out = double(temp_out)/255;
abs_out = exp(img_out*m) - 1;
fourier_out = abs_out.*exp(1i.*an);
img_prepared = ifft2(ifftshift(fourier_out));
%imwrite(img_prepared, "edited.png");
imshow(img_prepared, []);