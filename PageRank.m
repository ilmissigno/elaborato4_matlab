function [R,OUT,IN] = PageRank(G)
%% Algoritmo di PageRank
% L'algoritmo di Pagerank permette di ricavare il fattore di popolarita' di
% una pagina web (PageRank). Questo algoritmo data una matrice di
% adiacenza dei nodi di collegamento in un sottoinsieme di siti web G, nella forma di una
% matrice sparsa quadrata di valori logical, restituisce il valore di
% OutDegree e InDegree che indicano il numero di link che puntano alla
% pagina i-esima e che escono dalla pagina i-esima rispettivamente.

%Controllo se i parametri di input sono consistenti con la definizione
if(nargin < 1)
    error('Err:LESS_INPUT_ARGS','Inserire la matrice G di Adiacenza');
end

n = size(G,1); %dimensione della matrice G

controllo_GMatrix(G); % Effettuo i controlli su G

% Parametri di ingresso costanti
TOL = 10^-7;
NMAX = 200;
p = 0.85;

G = findSelfLoops(G); %funzione che risolve il problema degli autoloop, pone la diagonale di G nulla.

% Ricaviamo i valori necessari per l'applicazione della formula
e=ones(n,1); %vettore unitario
R = ones(n,1)/n; %soluzione di partenza indicata come rank fittizio pari a 1/n
N = full(sum(G)); %vettore degli outdegree

z = ((1-p)/n)*(N~=0) + (N==0)/n; %per risolvere il problema dei dangling node e rank sink

N(N~=0)=1./N(N~=0); % Calcolo la matrice diagonale (la memorizzo come un vettore!)

toll = checkUnderFlowTOL(TOL,R); % Controllo la convergenza della soluzione R sfruttando il valore di
% tolleranza.

numiter = 0; % Contatore numero iterazioni
errore = Inf; % Errore assoluto di riferimento
pGD = p*G.*N; %costante, non viene ricalcolata nel ciclo
%% Criterio di arresto
while(errore>toll || numiter == NMAX)
    numiter = numiter+1;
    r_temp = R; %salvo il valore precedente per utilizzarlo nella formula
    
    %% RANK
    R = (pGD*r_temp) + (e*z*r_temp); %x^k+1 = (pGD)*x^k + e*(z*x^k)
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