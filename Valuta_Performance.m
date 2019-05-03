function Valuta_Performance()
%% Script dei Test di Correttezza
%In questo script valutiamo la correttezza dell'algoritmo di pagerank
%andando a confrontare il risultato ottenuto rispetto alla funzione
%centrality del MATLAB
%{
N.B : Nel seguente codice, la funzione centrality richiede in ingresso: 
Il grafo ricavato dalla matrice G di adiacenza omettendo gli AutoLoop presenti
Il valore costante di probabilità P utilizzato nell'algoritmo di PageRank
Il fattore di Tolleranza adoperato in fase di calcolo dei Rank.
%}
[U,G] = surfer('http://www.unina.it',100);
grafo = digraph(G','Omitselfloops');
performance = centrality(grafo,'pagerank','FollowProbability',0.85,'Tolerance',10^-7);
colonne = grafo.outdegree;
righe = grafo.indegree;
[performance,indice] = sort(performance,'descend');
id = indice(1:15);
outdegrees = colonne(indice(1:15));
indegrees = righe(indice(1:15));
siti_web = U(indice(1:15));
table(id,performance(1:15),outdegrees,indegrees,'RowNames',siti_web,'VariableNames',{'ID','rank','outdegree','indegree'})
end