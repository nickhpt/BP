function [slowtimeBins] = getSlowtimebins(j,v,t)
%Inputs:
%j: Output image azimuth indicies
%v: Along-track scene velocity
%t: Along-track time vector

N = length(j);
slowtimeBins = zeros(1,N);
for i = 1:N
    [~, bin] = min(abs((j(i))*1.5/v - t ));
    slowtimeBins(i) = bin;      
end

end

