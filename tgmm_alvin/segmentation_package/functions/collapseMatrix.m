function V = collapseMatrix(M)
% This function collapses a 2D matrix M into a column vector V, where the
% value of row r of V is equal to the mean of the nonzero elements of
% row r of M. Equivilant to taking the profile of a line whose
% width equals the width of the ray. Written by NSJ on 13172025.
    h = height(M);
    V = zeros(h, 1);
    for row = 1:h
        value = mean(nonzeros(M(row, :)));
        V(row, 1) = value;
    end
end