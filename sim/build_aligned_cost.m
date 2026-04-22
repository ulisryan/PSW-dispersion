function aligned_cost = build_aligned_cost(particle, linearly, MAX)
% BUILD_ALIGNED_COST  Build cos2theta lookup matrix using vectorized indexing.
%
% Replaces the original sequential for loop:
%   for index = 1:MAX
%       aligned_cost(index,:) = linearly(particle(index,10)+1, particle(index,11)+1, :);
%   end
%
% Strategy: reshape the 3-D linearly table into a 2-D matrix (nj*nm × nI),
% then use linear indexing to extract all MAX rows in a single operation.
%
% Inputs:
%   particle  — MAX×11 particle array (columns 10,11 hold j, |m| quantum numbers)
%   linearly  — nj × nm × nI lookup table loaded from linearly.mat
%   MAX       — number of particles
%
% Output:
%   aligned_cost — MAX × nI matrix of cos2theta values per particle

    [nj, nm, ~] = size(linearly);

    % Collapse (j, m) dimensions into a single linear dimension
    linearly_2d = reshape(linearly, nj * nm, []);  % (nj*nm) × nI

    % 1-indexed row/column into the (j, m) plane
    j_idx = particle(:,10) + 1;   % MAX×1
    m_idx = particle(:,11) + 1;   % MAX×1

    % Single linear index for each (j, m) pair
    lin_idx = sub2ind([nj, nm], j_idx, m_idx);  % MAX×1

    % Extract all rows at once — no loop required
    aligned_cost = linearly_2d(lin_idx, :);  % MAX × nI

end
