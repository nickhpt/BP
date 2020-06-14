function [yint] = sincInterp(x,t,cell,M)

%Inputs
%x:Input data
%t:Sample points for input data
%Cell:The range that should be interpolated to
%M:Length of sinc-kernel (one-sided), must be an even number

%Interpolator output sample points
Rbin = cell;
nT = [floor(Rbin)-(M-1):floor(Rbin-1) floor(Rbin) Rbin ceil(Rbin) ceil(Rbin)+1:ceil(Rbin)+(M-1)];

yp = (zeros(M*2+1,M*2));
for j=1:length(x)       
    yp(:,j) = x(j)*sinc(t(j)-nT);
end

%Sum along second dimension to obtain interpolated data
yint = sum(yp,2);
%Return interpolated value
yint = yint(M+1);
end

