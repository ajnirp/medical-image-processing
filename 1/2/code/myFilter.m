%% Question 2 - Generating Filters
% if filterType = 1, means Ram-Lak filter
%
% if filterType = 2, means Shepp-Logan filter
%
% if filterType = 3, means Low pass cosine filter
%
% Length is as a factor of maximum frequency. So length of 0.75 means that
% all frequencies above 0.75 of maximum are assigned to zero.
function filterMatrix = myFilter(filterType, length)
if(nargin == 0)
    filterType = 3;
    length = 1;
end

if(nargin == 1)
    length = 1;
end

%%
% Ram Lak filter
if filterType == 1
    freqs = linspace(-1,1,367).';
    filterMatrix = abs(freqs);
    filterMatrix(abs(freqs)>length) = 0;
    filterMatrix = repmat(filterMatrix, [1 60]);
%     figure(); plot(filterMatrix); title('Ram Lak Filter');
end

%%
% Shepp Logan Filter
if filterType == 2
    freqs = linspace(-1,1,367).';
    filterMatrix = abs(freqs).*sin(pi*0.5*freqs)./(pi*0.5*freqs);
    filterMatrix(abs(freqs)>length) = 0;
    filterMatrix = repmat(filterMatrix, [1 60]);
%     figure(); plot(filterMatrix); title('Shepp-Logan Filter');
end

%%
% Low pass cosine filter
if filterType == 3
    freqs = linspace(-1,1,367).';
    filterMatrix = abs(freqs).*cos(pi*0.5*freqs);
    filterMatrix(abs(freqs)>length) = 0;
    filterMatrix = repmat(filterMatrix, [1 60]);
%     figure(); plot(filterMatrix); title('Low pass cosine filter')
end
