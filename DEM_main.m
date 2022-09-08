clc; clear all; close all;

% =========================================================================
% registers moving tissue map to fixed tissue map and outputs the computed 
% deformation field and warped tissue map
%
% Developer: Sahar Ahmad, 2022
% =========================================================================

% ------------------------------- %
%         path to data            %
% ------------------------------- %
directoryMoving = strcat('moving'); 
directoryFixed = strcat('fixed');
directoryOutput = strcat('output');


% ----------------------- %
%        subj ids         %
% ----------------------- %
MovingSubj = "Moving";
FixedSubj = "Fixed";

% ----------------------- %
%     initialization      %
% ----------------------- %
var = 0.4;

% ----------------------- %
%       open images       %
% ----------------------- %
filename = fullfile(directoryFixed, char(strcat(FixedSubj,'_tissue.nii.gz')));
if exist(filename, 'file') == 0
    disp('fixed map not found')
else
    info = niftiinfo(filename);
    I = double(niftiread(filename));
    size_x = size(I,1); size_y = size(I,2); size_z = size(I,3);
    
    I(I == 1) = 10;
    I(I == 2) = 150;
    I(I == 3) = 250;

    I = (I/max(I(:)));
    Ig = imgaussian(I, 3, var);
end

filename = fullfile(directoryMoving, char(strcat(MovingSubj,'_tissue.nii.gz')));
if exist(filename, 'file') == 0
    disp('moving map not found')
else
    temp = double(niftiread(filename));
    temp(temp == 1) = 10;
    temp(temp == 2) = 150;
    temp(temp == 3) = 250;

    I2 = temp;
    I2 = I2/(max(I2(:)));
    I2g = imgaussian(I2, 3, var);
    
    I22 = temp;
    mu = 1.1*ones(size(I22));
    mu(I22 == 10) = 1.1;
    mu(I22 == 150) = 0.9;
    mu(I22 == 250) = 0.8;
end


% ------------------------ %
ymin = 1; xmin = 1; zmin = 1;
ymax = size_x; xmax = size_y; zmax = size_z;

h = size_x;
w = size_y;
n = size_z;
img_res = 1.0;
[yt, xt, zt] = meshgrid(double(0:img_res:img_res*w - img_res), double(0:img_res:img_res*h - img_res), double(0:img_res:img_res*n - img_res));

Tinv = inv(info.Transform.T);
x = info.Transform.T(1,1).*xt;
y = info.Transform.T(2,2).*yt + info.Transform.T(3,2).*zt;
z = info.Transform.T(2,3).*yt + info.Transform.T(3,3).*zt;

ct = 0.001;

k = 5.5*ones(size(I));
rho = 5;
lambda = 0;

% ---------------- initial dice ratio --------------- %
S = 250*I2;
G = I*250;
listLabelS = unique(S);             % a list of labels of objects in S
listLabelS(listLabelS == 0) = [];

listLabelG = unique(G);             % a list of labels of objects in S
listLabelG(listLabelG == 0) = [];

for iLabelS = 1:length(listLabelS)
    Si = S == listLabelS(iLabelS);
    Gi = G == listLabelG(iLabelS);
    
    D(iLabelS) = 2*nnz(Si & Gi)/(nnz(Si) + nnz(Gi));
end
fprintf("Initial Dice Ratio: ")
mean(D)

lap = zeros(3, 3, 3);
lap(:, :, 1) = [0 1 0;1 2 1;0 1 0];
lap(:, :, 2) = [1 2 1;2 -24 2;1 2 1];
lap(:, :, 3) = [0 1 0;1 2 1;0 1 0];
lap = (1/6) * lap;

waterLevelPrevx = zeros(h, w, n);
waterLevelPrevy = zeros(h, w, n);
waterLevelPrevz = zeros(h, w, n);
U = zeros(h, w, n);
V = zeros(h, w, n);
W = zeros(h, w, n);

h1x = zeros(3, 3, 3);
h2x = h1x; h3x = h1x;
h1y = h1x; h2y = h1x; h3y = h1x;
h1z = h1x; h2z = h1x; h3z = h1x;
h1x(:, :, 2) = [0 1 0;0 -2 0;0 1 0];
h2x(:, :, 2) = 0.25 * [1 0 -1;0 0 0;-1 0 1];
h3x(:, :, 1) = 0.25 * [0 1 0;0 0 0;0 -1 0];
h3x(:, :, 3) = 0.25 * [0 -1 0;0 0 0;0 1 0];
h1y(:, :, 2) = 0.25 * [1 0 -1;0 0 0;-1 0 1];
h2y(:, :, 2) = [0 0 0;1 -2 1;0 0 0];
h3y(:, :, 1) = 0.25 * [0 0 0;1 0 -1;0 0 0];
h3y(:, :, 3) = 0.25 * [0 0 0;-1 0 1;0 0 0];
h1z(:, :, 1) = .25 * [0 1 0;0 0 0;0 -1 0];
h1z(:, :, 3) = .25 * [0 -1 0;0 0 0;0 1 0];
h2z(:, :, 1) = .25 * [0 0 0;1 0 -1;0 0 0];
h2z(:, :, 3) = .25 * [0 0 0;-1 0 1;0 0 0];
h3z(:, :, 1) = [0 0 0;0 1 0;0 0 0];
h3z(:, :, 2) = [0 0 0;0 -2 0; 0 0 0];
h3z(:, :, 3) = [0 0 0;0 1 0;0 0 0];

f = zeros([size(I,1), size(I,2), size(I,3), 3]);
step = 1e-05;
trm = 30; 
Im = I;

for loop = 1:21
    cnt = 0; chk = 0;
    fprintf("Level: ")
    loop
    for i = 1:trm
        X = xt + U;
        Y = yt + V;
        Z = zt + W;
        
        [I2st, ~] = interpolation3D(I2g, xt, yt, zt, X, Y, Z);
        
        graddivx = imfilter(U, h1x, 'replicate') + imfilter(V, h2x, 'replicate') + imfilter(W, h3x, 'replicate');
        graddivy = imfilter(U, h1y, 'replicate') + imfilter(V, h2y, 'replicate') + imfilter(W, h3y, 'replicate');
        graddivz = imfilter(U, h1z, 'replicate') + imfilter(V, h2z, 'replicate') + imfilter(W, h3z, 'replicate');
        lapimg1 = imfilter(U, lap, 'replicate');
        lapimg2 = imfilter(V, lap, 'replicate');
        lapimg3 = imfilter(W, lap, 'replicate');
        
        
        [I2gx, I2gy, I2gz] = gradient(I2st);
        
        force = Im - I2st;
        fx = k .* force .* (I2gy);
        fx(isnan(fx)) = 0;
        
        fy = k .* (force) .*(I2gx);
        fy(isnan(fy)) = 0;
        
        fz = k .* (force) .*(I2gz);
        fz(isnan(fz)) = 0;
        
                
        fx = imgaussian(fx, 3, 0.8);
        fy = imgaussian(fy, 3, 0.8);
        fz = imgaussian(fz, 3, 0.8);

        
        if (abs(mean(fx(:))) > step)
            waterLevelNextx = ( ( mu .* lapimg1  + ( mu + lambda ) .* graddivx + fx)/rho ) + 2 * U - waterLevelPrevx;
            flagx = 0;
        else
            waterLevelNextx = U;
            flagx = 1;
        end
        
        if (abs(mean(fy(:))) > step)
            waterLevelNexty = ( ( mu .* lapimg2  + ( mu + lambda ) .* graddivy + fy)/rho ) + 2 * V - waterLevelPrevy;
            flagy = 0;
        else
            waterLevelNexty = V;
            flagy = 1;
        end
        
        if (abs(mean(fz(:))) > step)
            waterLevelNextz = ( ( mu .* lapimg3  + ( mu + lambda ) .* graddivz + fz)/rho ) + 2 * W - waterLevelPrevz;
            flagz = 0;
        else
            waterLevelNextz = W;
            flagz = 1;
        end
        
        if flagx == 1 && flagy == 1 && flagz == 1
            break
        end
        
        waterLevelPrevx = U;
        U = waterLevelNextx;
        
        waterLevelPrevy = V;
        V = waterLevelNexty;
        
        waterLevelPrevz = W;
        W = waterLevelNextz;
        k = k + .08;
    end
    
    ct = ct + 0.002;
    if mu > 0.3
        indt = find(I22 == 150);
        mu(indt) = mu(indt) - 0.01;
        indt = find(I22 == 250);
        mu(indt) = mu(indt) - 0.01;
    end
    
    U = waterLevelPrevx;
    V = waterLevelPrevy;
    W = waterLevelPrevz;
    
    X = xt + U;
    Y = yt + V;
    Z = zt + W;
    
    [I2st, ~] = interpolation3D_seg(I2, xt, yt, zt, X, Y, Z);
    S = 250*I2st;
    G = I*250;
    listLabelS = unique(S);             % a list of labels of objects in S
    listLabelS(listLabelS == 0) = [];
    
    listLabelG = unique(G);             % a list of labels of objects in S
    listLabelG(listLabelG == 0) = [];
    
    for iLabelS = 1:length(listLabelS)
        Si = S == listLabelS(iLabelS);
        Gi = G == listLabelG(iLabelS);
        
        D(iLabelS) = 2*nnz(Si & Gi)/(nnz(Si) + nnz(Gi));
    end
    
    fprintf("Dice Ratio is ")
    mean(D)
    
end

% -------------- Save output -------------------- %
deformation = zeros(size_x,size_y,size_z,1,3);
for i = 1:size(U,3)
deformation(ymin:ymax,xmin:xmax,zmin + i - 1,1,1) = U(:,:,i);
deformation(ymin:ymax,xmin:xmax,zmin + i - 1,1,2) = V(:,:,i);
deformation(ymin:ymax,xmin:xmax,zmin + i - 1,1,3) = W(:,:,i);
end

niftiwrite(deformation,fullfile(directoryOutput,char(strcat('deform_',MovingSubj,'_To_',FixedSubj,'_DEMWarp.nii'))));

[I2st, ~] = interpolation3D_seg(I2, xt, yt, zt, X, Y, Z);

I1 = zeros(size_x, size_y, size_z);
I1(ymin:ymax, xmin:xmax, zmin:zmax) = I2st;
Im = zeros([size(I1,1),size(I1,2),size(I1,3)]);
arr = zeros(3,3,3);
for k = 2:size(I1,3)-1
    for j = 2:size(I1,2)-1
        for i = 2:size(I1,1)-1
            arr = I1(i - 1:i + 1,j - 1:j + 1,k - 1:k + 1);
            vote = mode(arr(:));
            Im(i,j,k) = vote;
        end
    end
end
    
filename = fullfile(directoryOutput, char(strcat(MovingSubj,'_To_',FixedSubj,'_tissue_reg_DEM.nii')));
niftiwrite(250*Im,filename,info);

