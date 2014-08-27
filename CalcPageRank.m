function [val,rank] = CalcPageRank(mat);

Centrality = pagerank(mat);
Centrality = (Centrality - min(Centrality))./(max(Centrality)-min(Centrality));
[val,rank]=sort(Centrality,'descend');

    function x = pagerank(G)
        % PAGERANK  Google's PageRank
        % pagerank(U,G,p) uses the URLs and adjacency matrix produced by SURFER,
        % together with a damping factory p, (default is .85), to compute and plot
        % a bar graph of page rank, and print the dominant URLs in page rank order.
        p = .85;
        % Eliminate any self-referential links
        G = G - diag(diag(G));
        % c = out-degree, r = in-degree
        [~,n] = size(G);
        c = full(sum(G,1));
        r = full(sum(G,2));
        % Scale column sums to be 1 (or 0 where there are no out links).
        k = find(c~=0);
        D = sparse(k,k,1./c(k),n,n);
        % Solve (I - p*G*D)*x = e
        e = ones(n,1);
        I = speye(n,n);
        x = (I - p*G*D)\e;
        % Normalize so that sum(x) == 1.
        x = x/sum(x);
    end
end