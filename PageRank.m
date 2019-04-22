function [R,OUT,IN] = PageRank(G)
%% Algoritmo di PageRank
% L'algoritmo di Pagerank permette di ricavare il fattore di popolarita' di
% una pagina web (PageRank). Questo algoritmo data una matrice di
% adiacenza dei nodi di collegamento in un sottoinsieme di siti web G, nella forma di una
% matrice sparsa quadrata di valori logical, restituisce il valore di
% OutDegree e InDegree che indicano il numero di link che puntano alla
% pagina i-esima e che escono dalla pagina i-esima rispettivamente.

% Per prima cosa bisogna ottenere la dimensione della matrice per risolvere
% i problemi di autoloop, rank-sink e dangling node.
% Il problema di autoloop e' relativo a pagine che si autoreferenziano
% ovvero producono outdegree che puntano a se stessi.
% Il problema del rank-sink e' relativo a pagine che hanno un solo link
% entrante e nessun link uscente (es. problema della trappola) ed hanno un
% rank sempre piu' grande dato che il grafo non e' fortemente connesso.
% Il problema dei dangling node e' relativo al problema del rank-sink ma le
% pagine hanno rank uguale a 0.

%Controllo se i parametri di input sono consistenti con la definizione
if(nargin < 1)
    error('Err:LESS_INPUT_ARGS','Inserire la matrice G di Adiacenza');
end

% Ricavo la dimensione
n = size(G,1);
% Effettuo i controlli di routine su G
controllo_GMatrix(G);

% Per determinare il risultato in maniera corretta, l'algoritmo deve
% lavorare su un criterio di arresto. Quest'ultimo deve essere imposto
% andando a verificare che l'errore di round-off relativo ai dati di input
% non deve superare un valore di Tolleranza relativo all'ultima soluzione
% trovata in fase di calcolo del rank (sfruttando il metodo iterativo
% imposto dall'algoritmo di pagerank).
% La formula che il ciclo, veicolato dal criterio di arresto e dalla
% condizione di emergenza, che utilizza per calcolare il rank e' la
% seguente.
%% Rank : x^k+1 = (p*G.*N)*x^k + e*(z*x^k)
% Questa formula fa uso del valore all'istante k del rank moltiplicato il
% vettore z colonna indica i possibili dangling nodes che sono presenti nella
% matrice G, il vettore 'e' e il fattore correttivo utilizzato per evitare
% presenze di rank-sink, 'p=0.85' e' la costante moltiplicativa utilizzata
% per dividere parte del rank sui link uscenti, ed N e' il vettore di
% outdegree ottenuto sommando tutti gli elementi di G e costruendone un
% vettore pieno (ond'evitare puntamenti a pagine non presenti).

% Parametri di ingresso costanti
TOL = 10^-7;
NMAX = 200;
p = 0.85;

% Prima di ricavare i valori necessari per il calcolo del pagerank, bisogna
% risolvere il problema degli autoloop, ponendo la diagonale di G nulla per
% eliminare link a se stessi.

G = findSelfLoops(G);

% Ricaviamo i vari valori necessari per l'applicazione della formula

e=ones(n,1); %deve essere un vettore unitario ond'evitare rank-sink per ogni nodo.

%prima di ricavare gli altri parametri e' necessario utilizzare una
%soluzione di partenza indicata come rank fittizio pari a 1/n
R = ones(n,1)/n;

%calcoliamo il vettore degli outdegree
N = full(sum(G));

%Per ricavare il vettore z bisogna tenere conto del problema del rank-sink.
%Ovviamente per uscire dalla situazione del rank sink bisogna sceglie a
%caso il valore di un link uscente oppure saltare ad altre pagine
%arbitrarie. dunque la probabilita' che si effettui la prima scelta viene
%indicata con p, di conseguenza l'altra e' indicata con probabilita' 1-p.
%In definitiva, se p/Nj e' la probabilita' di scegliere un link uscente su
%una pagina j, e (1-p)/Nj e' la probabilita' di saltare ad una pagina j
%arbitraria, allora ogni valore di z puo' essere calcolato come:
%% zj = p(1/Nj) + (1-p)/n per Nj ~= 0
%% zj = 1/n per Nj == 0 (dangling node salto alla pagina j arbitraria)

%uso la notazione di matlab che automaticamente scorre il vettore N degli
%outdegree e controlla se ci sono elementi nulli o meno per poi
%moltiplicarlo al valore di probabilita' corrispettivo e ricava un vettore
%colonna.
z = ((1-p)/n)*(N~=0) + (N==0)/n;

% Calcolo la matrice diagonale (la memorizzo come un vettore!)
N(N~=0)=1./N(N~=0);

% Controllo la convergenza della soluzione R sfruttando il valore di
% tolleranza.
toll = checkUnderFlowTOL(TOL,R);

numiter = 0; % Contatore numero iterazioni
errore = Inf; % Errore assoluto di riferimento

%% Criterio di arresto
while(errore>toll || numiter == NMAX)
    numiter = numiter+1;
    r_temp = R; %salvo il valore precedente per utilizzarlo nella formula
    
    %% RANK
    R = ((p*G.*N)*r_temp) + (e*z*r_temp);
    %ricalcolo l'errore (assoluto) sulla soluzione attuale
    errore = norm(R-r_temp,Inf);
    %ricalcolo la convergenza
    toll = checkUnderFlowTOL(TOL,R);
end

%% Ricavo degli InDegree e OutDegree da G
IN = full(sum(G,2));
OUT = full(sum(G,1))';

end

function controllo_GMatrix(G)
%% Funzione di controllo della Matrice G
% La matrice G deve essere una matrice sparsa, quadrata e logica nel senso
% che deve avere valori booleani o semplicemente 0 e 1.

if(~ismatrix(G))
       error('Err:G_isMatrix','G deve essere una matrice');
end

if(isempty(G))
       error('Err:EmptyMatrix','La matrice G non deve essere vuota');
end
if(find(~isfinite(G)))
    error('Err:Matrix_NotNumeric','I valori di G devono essere numeri reali finiti');
end
if(size(G,1)~=size(G,2))
     error('Err:Matrix_NotQuadratic','La matrice G di ingresso deve essere quadrata');
end
if(~issparse(G))
      error('Err:Matrix_NotSparse','La matrice G deve essere sparsa');
end

if(~islogical(G))
    error('Err:G_NON_LOGICAL','La matrice G deve avere soltanto valori booleani o valori pari a 0 o 1');
end
end

function G = findSelfLoops(G)
%% Funzione di eliminazione degli autoloop della Matrice G
if(any(find(abs(diag(G)) ~= 0)))
    warning('Warn:SELF_LOOPS_EXISTS','Presenti selfloops, ne effettuo la rimozione');
    G = spdiags(zeros(size(G,1),1),0,G);
end
end

%% Funzione di verifica del Teorema di Convergenza del metodo iterativo
function toll = checkUnderFlowTOL(TOL,solk)
    toll = TOL*norm(solk,Inf);
    if (toll < realmin)
        toll = realmin;
    end
end