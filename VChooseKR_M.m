function Y = VChooseKR_M(X, K)
% VChooseKR - Combinations of K elements, repetitions [Matlab]
% See VChooseKR for help.
% In Matlab 7.8 the function run fastest for DOUBLEs.
%
% Surprising: VChooseK_M is faster than COMBINATOR for small K and large N, e.g.
% (int16(1:2048), 2): 35 times faster than COMBINATOR!
% But in general the speed of COMBINATOR exceeds this Matlab-version remarkable,
% e.g. (uint8(1:2), 64): COMBINATOR is 32 times faster than VChooseK_M.
% Of course, VChooseK.mex should be used for speed in real calculations.
%
% Tested: Matlab 6.5, 7.7, 7.8, WinXP, [UnitTest]
% Author: Jan Simon, Heidelberg, (C) 2009-2010 matlab.THISYEAR(a)nMINUSsimonDOTde
% License: BSD (use/copy/modify on own risk, but mention author)

% ==============================================================================
% Slower M-version as proof of concept - no input checks here.
X  = X(:);         % Give X a column shape
nX = numel(X);
Y  = X([]);        % Duplicate class of the input
if nX == 0 || K == 0
    return;
end

switch K
    case 1
        Y = X;
        
    case 2
        % Matlab 7: Y = zeros(nY, 3, class(X));
        nY       = round(((nX + 1) * nX) / 2);
        Y(nY, 2) = 0;  % Class of input is kept
        a  = 1;
        for k = 1:nX
            b         = a + nX - k;
            Y(a:b, 1) = X(k);
            Y(a:b, 2) = X(k:nX);
            a         = b + 1;
        end
        
    case 3  % Just insert another loop:
        % Matlab 7: Y = zeros(nY, 3, class(X));
        nY       = round(((nX + 2) * (nX + 1) * nX) / 6);
        Y(nY, 3) = 0;
        a  = 1;
        for k1 = 1:nX
            for k2 = k1:nX
                b         = a + nX - k2;
                Y(a:b, 1) = X(k1);
                Y(a:b, 2) = X(k2);
                Y(a:b, 3) = X(k2:nX);
                a         = b + 1;
            end
        end
        
    case 4   % Insert another loop again:
        % 4.4 times faster than Matlab's NCHOOSEK, but the MEX is 240 times faster
        % for (1:20).
        % Matlab 7: Y = zeros(nY, 3, class(X));
        nY = round(((nX + 3) * (nX + 2) * (nX + 1) * nX) / 24);
        Y(nY, 4) = 0;
        a = 1;
        for k1 = 1:nX
            for k2 = k1:nX
                for k3 = k2:nX
                    b         = a + nX - k3;
                    Y(a:b, 1) = X(k1);
                    Y(a:b, 2) = X(k2);
                    Y(a:b, 3) = X(k3);
                    Y(a:b, 4) = X(k3:nX);
                    a         = b + 1;
                end
            end
        end
        
    otherwise
        % It is trivial to expand this more and more by inserting loops over
        % k4, k5, k6... The general algorithm is more challenging:
        
        % Calculate number of rows and create output:
        nY       = NoverK(nX + K - 1, K);
        Y(nY, K) = 0;  % Keep type of input
        
        % The algorithm is designed in C-style, because it was used to debug the
        % C-Mex function. Just the loop over the K.th column is vectorized. With
        % using the power of Matlab, a much faster algorithm would be possible
        % (see e.g. COMBINATOR of Matt Fig, FEX: 24325). But this is just a dummy
        % for the much faster MEX!
        Index    = ones(1, K);             % Current index
        Limit    = nX(ones(1, K));
        Limit(K) = 0;
        a        = 1;
        while 1
            b = a + nX - Index(K);          % Write index for last column
            for i = 1:(K - 1)               % Write the left K-1 columns
                Y(a:b, i) = X(Index(i));
            end
            Y(a:b, K) = X(Index(K):nX);     % Write the K.th column
            a         = b + 1;              % Move the write pointer
            
            % Search the last column index, which is not exhausted:
            newLoop = find(Index < Limit);
            if isempty(newLoop)             % All columns are filled:
                return;                      % Ready!
            end
            newLoop = newLoop(length(newLoop));
            
            % Fill new Index with new value encreasing by 1:
            Index(newLoop:K) = Index(newLoop) + 1;
        end
end  % switch K

return;

    function m = NoverK(n, k)
        % Use a loop to calculate the length of the output. All intermediate values are
        % integers to reduce rounding problems.
        
        % (n over k) == (n over n-k), so use the shorter sequence:
        if k + k > n
            k = n - k;
        end
        
        % Catch K==0 and K==1:
        if k <= 1
            if k == 1
                m = n;
            else
                m = 1;
            end
            return;
        end
        
        % Idea:
        % Intermediate values are not integers:
        %   m = prod(((n-k+1):n) ./ (1:k))
        % Intemediate values are integers:
        %   m = prod((...((n-k+1) * (n-k+2) / 2) * (n-k+3) / 3) * ... * n / k)
        % Calculate (n-k+1) * (n-k+2) / 2:
        %   m_1 = (n-k+1)
        %   m_2 = (n-k+1) * (n-k+2) / 2 =
        %       = (n-k+1) + (n-k+1) * (n - k) / 2 = m_1 + m_1 * (n - k) / 2
        d = n - k;
        m = d + 1;
        for i = 2:k
            m = m + (m * d) / i;
        end
        
        return;
    end
end
