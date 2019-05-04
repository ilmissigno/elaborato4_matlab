classdef CasiTest < matlab.unittest.TestCase
    
    methods(Test)
        
        function TestCase1(testCase)
            %% TEST NUMERO PARAMETRI IN INGRESSO
            %verifica l'errore nel caso in cui il numero di parametri di
            %ingresso e' minore di 1
            verifyError(testCase,@()PageRank(),'Err:LESS_INPUT_ARGS');
        end
        
        function TestCase2(testCase)
            %% TEST MATRICE G VUOTA
            % Verifica l'errore nel caso in cui la matrice G è vuota
            G = Richiama_Parametri();
            G = [];
            verifyError(testCase,@()PageRank(G),'Err:EmptyMatrix');
        end
        
        function TestCase3(testCase)
            %% TEST MATRICE G NON SPARSA
            %Verifica l'errore nel caso in cui la matrice G non è sparsa
            
            G = Richiama_Parametri();
            G = rand(10,10);
            verifyError(testCase,@()PageRank(G),'Err:Matrix_NotSparse');
            
        end
        
        function TestCase4(testCase)
            %% TEST MATRICE G NON FINITA
            %verifica l'errore nel caso in cui la matrice G non è finita
            G = Richiama_Parametri();
            G = [2 3; inf 4];
            verifyError(testCase,@()PageRank(G),'Err:Matrix_NotNumeric');
            
        end
        
        function TestCase5(testCase)
            %% TEST MATRICE G NON QUADRATICA
            %verifica l'errore nel caso in cui la matrice G non è
            %quadratica
            G = Richiama_Parametri();
            G = [2 3 4; 4 5 3];
            verifyError(testCase,@()PageRank(G),'Err:Matrix_NotQuadratic');
            
        end
        
        function TestCase6(testCase)
            %% TEST MATRICE G NON LOGICA
            %verifica l'errore nel caso in cui G non e' logical
            G = Richiama_Parametri();
            G = round(sprand(10,10,0.2)*10);
            verifyError(testCase,@()PageRank(G),'Err:G_NON_LOGICAL');
        end
        
        
        
        function TestCase7(testCase)
            %% TEST SELF LOOPS
            %verifica se sono presenti self loops all'interno della matrice
            G = Richiama_Parametri();
            G = logical(G+speye(size(G)));
            verifyWarning(testCase,@()PageRank(G),'Warn:SELF_LOOPS_EXISTS');
        end
    end
end

function G = Richiama_Parametri()
G = sparse(logical(randi([0,1],10,10)));
end