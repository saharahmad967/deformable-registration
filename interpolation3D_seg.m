function  [Interp_Image,Mask] = interpolation3D_seg(imgI,X,Y,Z,Fx,Fy,Fz)

Interp_Image = interp3(Y,X,Z,double(imgI),Fy,Fx,Fz,'nearest');
ind = find(isnan(Interp_Image));
Interp_Image(ind) = 0;
Mask = ones(size(Interp_Image));  
Mask(ind) = 0;
