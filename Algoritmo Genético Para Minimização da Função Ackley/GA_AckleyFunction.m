%---------------------------------------------------------------------------------------------------------------------------------------%
    % Trabalho 03 - Inteligência Computacional - Algoritmos Geneticos
    % Emanuel Valerio Pereira - Matricula 471055
%---------------------------------------------------------------------------------------------------------------------------------------%

% Intervalo de x e y onde a funcao foi definida
x = -10:0.1:10;
y = -10:0.1:10;

% Plotando a funcao em 3D, apenas para fins de vizualizacao
[X Y] = meshgrid(x,y);
z = (-20*exp(-0.2*sqrt(0.5*((X).^2+(Y).^2))) - exp(0.5*(cos(2*pi*X)+cos(2*pi*Y)))+exp(1) + 20);
 figure,surf(x,y,z),hold on
axis('tight'),xlabel('x'),ylabel('y'),title('Peaks')
randomPoints = randperm(numel(z));

ind = 50; % Numero de Individuos ( numero de cromossomos );

for i = 1:ind
    plot3(X(randomPoints(i)),Y(randomPoints(i)),z(randomPoints(i)),'r.','MarkerSize',20),
    hold on
end

% Definicao de alguns parametros importantes
maxGen = 350; % Numero maximo de geracoes
numVar = 2; % numero de variaveis, nesse problema x e y.
numBits = 20; % Numero de Genes ( Numero de Bits ).
txm = 0.05; % Taxa de Mutação 5%
gen = 1; % Geração inicial
numDeath = round( 0.2*ind); % Numero maximo de inviduos que nao avancaram para a proxima geracao

% Inicializacao de algumas matrizes, usadas na resolucao do problema
cromossomos = zeros(ind,numBits*numVar); 
fitness = zeros(ind,1);
nota = zeros(ind,1);
notaDeath = zeros(ind,1);
Real = zeros(ind,numVar);
RealDeath = zeros(ind,numVar);
newGen = zeros(ind,numBits*numVar);
notDeath = zeros(40,1);
prob = zeros(ind,1);
probDeath = zeros(ind,1);
prob_roleta = zeros(ind,1);
prob_roleta_death = zeros(ind,1);
prob_acumulada = zeros(ind,1);
prob_acumulada_death = zeros(ind,1);
cromossomosAux = zeros(ind,40);

% Gerando a populacao(cromossomos) inicial 
cromossomos = randi([0 1],ind,numVar*numBits);

while gen<=maxGen
% Inicializando as variaveis contadoras 
    newInd=1;
    nota = 0;
    sumFitness = 0;
    sumNota = 0;
% Calculando o valor real correspodente a cada individuo da populacao, uma
% vez que estao dispostos em binario, abaixo eh feito a conversao para
% numero real.
    for i = 1:ind
        Real(i,1) = min(x) + (bin2dec(int2str(cromossomos(i,1:20)))/(2^20-1))*(max(x)-min(x));
        Real(i,2) = min(y) + (bin2dec(int2str(cromossomos(i,21:40)))/(2^20-1))*(max(y)-min(y));
    end
% Uma vez calculado o valor real, abaixo esta sendo feita a avaliacao de
% cada individuo, (atribuicao de notas).
    for i = 1:ind
           nota(i)=(-20*exp(-0.2*sqrt(0.5*((Real(i,1)).^2+(Real(i,2)).^2))) - exp(0.5*(cos(2*pi*Real(i,1))+cos(2*pi*Real(i,2))))+exp(1) + 20); 
           sumNota = sumNota + nota(i); % Realizando a soma de todas as notas (avaliacao)
    end
% O problema trata-se de uma minimizacao, logo, a menor nota obtida no
% calculo anterior, sera o melhor individuo e deve possuir maior
% probabilidade de cruzamento, entao, se apenas dividir-mos a nota pela
% soma, ao contrario do que queremos, o melhor individuo tera a menor
% probabilidade, para corrigir as notas de tal forma de contornar esse
% problema, temos :

    for i = 1:ind
           fitness(i) = ((max(nota)-nota(i))+min(nota));
           sumFitness = sumFitness+fitness(i); % Soma de todas as notas corrigidas (fitness)
    end
% Com as avaliacoes corrigidas, pode-se entao calcular a probabilidade de
% cada individuo ser selecionado para cruzamento.
    for i = 1:ind
        prob(i) = fitness(i)/sumFitness;
    end
% Para simular a roleta, multiplica-se as probabilidades por 360, com esse
% resultado teremos a fatia da roleta (angulo) correspodente a cada
% individuo.
    prob_roleta = prob*360.0;

% Gerando a roleta virtual ( calculo da Probabilidade acumulada )
    indAux = 0;
    for i = 1:ind
        indAux = indAux + prob_roleta(i);
        prob_acumulada(i) = indAux;
    end
% Selecao dos pais (Individuos para cruzamento)
while newInd <= ind
   P = 1 + (360-1).*rand(1,1); % Gerando um numero randomico entre 1 e 360 (angulo roleta)
        for j = 1:ind
            if prob_acumulada(j) >= P 
                pai1 = j; % Selecionando o Pai1
                break;
            end
        end
        % O indice do individuo selecionado eh o primeiro valor da
        % probabilidade acumulada, que supere o valor aleatorio gerado
   P = 1 + (360-1).*rand(1,1);
       for j = 1:ind
            if prob_acumulada(j) >= P
                pai2 = j; % Selecionando o Pai2
                break;
            end
       end 

  % Realizando o Crossover
        c=  randi([1 (numBits*(numVar-1))],1,1); % Gera um numero aleatorio entre 1 e o tamanho do cromossomo-1.
        gene11 = cromossomos(pai1,1:c);
        gene12 = cromossomos(pai1,c+1:numVar*numBits);
        gene21 = cromossomos(pai2,1:c);
        gene22 = cromossomos(pai2,c+1:numVar*numBits);
        filho1 = [gene11 gene22]; % Concatenacao dos genes para gerar filho1
        filho2 = [gene21 gene12]; % Concatenacao dos genes para gerar filho2

   % Armazenando os novos individuos em newGen, representando a nova geracao
        newGen(newInd,:)=  filho1;
        newInd = newInd+1;
        newGen(newInd,:)=  filho2;
        newInd = newInd + 1;

  % Operador de Mutacao

        if txm > rand() % se a taxa de mutacao for superior a um numero aleatorio entre 0 e 1, entao
            d = round(1+(numBits*numVar-2)*rand(1,1)); % Gera um valor aleatorio correspodente a um Bit (Gene) do cromossomo
            if newGen(newInd-2,d) == 0 
                newGen(newInd-2,d) = 1;
            else
                newGen(newInd-2,d) = 0;
            end
        % Como newInd esta sendo incrementada 2 vezes, faz-se necessário
        % voltar os indices em -2 e -1 para acessar os filhos gerados nesse
        % laço, que eh o filho1 e filho2

            if newGen(newInd-1,d) == 0
                newGen(newInd-1,d) = 1;
            else
                newGen(newInd-1,d) = 0;
            end
        end
end

% Aplicacao do conceito de Eletismo

    fitnessAux = zeros(ind,1); % Inicializa uma matriz auxiliar
    notaAux = zeros(ind,1);    % Inicializa uma matriz de notas auxiliar
    for i = 1:ind
        fitnessAux(i) = fitness(i); % Copia os valores de fitness para a matriz auxiliar
        notaAux(i) = nota(i);       % Copia as avaliacoes para a matriz auxiliar 
    end

  deathInd = zeros(numDeath,1); % Inicializa a matriz que ira armazenar os indices dos individuos com piores notas da novaGeracao
  cromossomosAux = cromossomos; % Copia dos individuos em uma matriz auxiliar
  sumNota=0;
  % Como queremos encontrar os piores individuos gerados para os eliminar,
  % e dessa forma resgatar a mesma quantidade de individuos so que os
  % melhores da geracao passada, para isso sera necessario criar uma roleta
  % da morte, dessa vez a probabilidade se da pela funcao inicial (nota), pois ela
  % ira penalizar os piores individuos ja que o objetivo agora eh eliminar
  % individuos, e aqueles com maiores notas, terao maior probabilidade de
  % serem eliminados.

  for i = 1:ind
         RealDeath(i,1) = -10 + (bin2dec(int2str(newGen(i,1:numBits)))/((2^20)-1))*(10-(-10));
         RealDeath(i,2) = -10 + (bin2dec(int2str(newGen(i,numBits+1:2*numBits)))/((2^20)-1))*(10-(-10));
   end
    for i = 1:ind
        notaDeath(i)=(-20*exp(-0.2.*sqrt(0.5*(((RealDeath(i,1)).^2)+((RealDeath(i,2)).^2)))) - exp(0.5*(cos(2.*pi.*Real(i,1))+cos(2.*pi.*Real(i,2))))+exp(1) + 20); 
        sumNota = sumNota+notaDeath(i);  
    end
    % Calculo da probabilidade de cada individuo gerado ser eliminado
    for i = 1:ind
        probDeath(i) = notaDeath(i)/sumNota;
    end
    % Multiplicando por 360, correspodente as fatias (angulo) que cada
    % individuo ira ocupar.
    prob_roleta_death = probDeath*360.0;

  % gerando a roleta ( calculo da Probabilidade acumulada) 
    indAux = 0;
    for i = 1:ind
        indAux = indAux + prob_roleta_death(i);
        prob_acumulada_death(i) = indAux;
    end
 % Selecao dos individuos que serao eliminados;
 % O indice do individuo selecionado eh o primeiro valor da
 % probabilidade acumulada, que supere o valor aleatorio (P) gerado.
    for i = 1:numDeath
     P = 1 + (360-1).*rand(1,1);
        for j = 1:ind
            if prob_acumulada_death(j) >= P
                indDeath = j;
                break;
             end
        end
        deathInd(i) = indDeath; % Armazenando os indices dos individuos que serao eliminados
    end
 % Como a selecao pode ocorrer com individuos repetidos, estou usando a
 % funcao unique para eliminar repeticoes, e assim numDeath sera um valor
 % maximo de individuos selecionados (caso nao haja repeticoes) caso
 % exista, esse valor de individuos a serem eliminados sera menor.
    deathInd = unique(deathInd);
    [li cl] = size(deathInd);
     best = zeros(li,40); % Inicializando a matriz que armazenara os melhores individuos da geracao passada
  
 % Eliminando os individuos selecionados 
    for i = 1:length(deathInd)
      newGen(deathInd(i),:) = [];
      deathInd = deathInd-1;
    end
 % A nova geracao possui agora length(newInd) a menos, uma vez que foram
 % eliminados, abaixo eh feita a atualizacao dos cromossomos com a nova
 % geracao.

    cromossomos = newGen;
 % Resgatar os melhores individuos da geracao anterior, nesse caso como
 % esses melhores individuos seram armazenados junto da nova geracao, o
 % numero de melhores individuos da antiga geracao deve ser igual ao
 % numero de individuos eliminados da nova.

      for b = 1:li
        [maxi indMax] = max(fitnessAux);
        best(b,:) = cromossomosAux(indMax,1:40); % Best armazena os 'li' melhores individuos
        fitnessAux(indMax) = 0;
      end

      cromossomos(ind-li+1:ind,:) = best(1:li,:); % Atualizando os individuos, juntando os individuos da nova geracao com os 
      % melhores individuos da antiga, garantindo que a proxima geracao
      % seja melhor que a atual
  
    [maior indice_maior] = max(fitness); % Indice do melhor individuo
    [menor indice_menor] = min(fitness); % indice do Pior individuo
% Plotagem atualizada do progresso do algoritmo

    contourf(x,y,z,'LineStyle','none'),hold on
    axis('tight'),xlabel('x','FontSize',16),ylabel('y','FontSize',16),
    title(sprintf('Geração %d/%d',gen,maxGen),"FontSize",20);
    for j = 1:length(Real)
    plot3(Real(j,1),Real(j,2),fitness(j),'r.','MarkerSize',20);
    end
    hold off
    pause(0.1)

    fprintf('============================================================\n');
    fprintf('                        GERACAO %d                          \n',gen);
    fprintf(' Nota Média : %.8f\n',sumFitness/ind);
    fprintf(' Nota Pior Individuo : %.8f , Porcetagem da Roleta : %.5f%% \n',nota(indice_menor),(100.*prob(indice_menor)));
    fprintf(' Nota Melhor Individuo : %.8f , Porcetagem da Roleta : %.5f%% \n',nota(indice_maior),(100.*prob(indice_maior)));
    if gen == maxGen
            fprintf('============================================================\n');
            fprintf('                        SOLUÇÃO                               \n');
            fprintf('Valor de X : %.6f\n',Real(indice_maior,1));
            fprintf('Valor de Y : %.6f\n',Real(indice_maior,2));

            figure,surf(x,y,z),hold on
            axis('tight'),xlabel('x'),ylabel('y'),title('Peaks')

            for i = 1:ind
                plot3(Real(i,1),Real(i,2),fitness(i),'r.','MarkerSize',20),
                hold on
            end
    end
    gen = gen +1 ;
end

% Menor nota obtida = 0.000038