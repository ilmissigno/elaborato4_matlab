[U,G] = surfer('http://www.unina.it',100);
grafo = digraph(G');
[R,OUT,IN] = PageRank(G);
grafo.Nodes.PageRank = R;
grafo.Nodes.InDegree = IN;
grafo.Nodes.OutDegree = OUT;
plot(grafo,'NodeColor','b','EdgeColor',[.3 .3 .3]);
axis off
%aaaa
