function f = imgaussian(f, siz, sigma)
x = -ceil(siz/2):ceil(siz/2);
% x = -0.5:0.1:0.5;
H = exp(-(x.^2/(2*sigma.^2)));
H = H/sum(H(:));
Hx = reshape(H, [length(H) 1 1]);
Hy = reshape(H, [1 length(H) 1]);
Hz = reshape(H, [1 1 length(H)]);
f = imfilter(imfilter(imfilter(f, Hx, 'same' ,'replicate'), Hy, 'same' ,'replicate'),Hz, 'same' ,'replicate');
end