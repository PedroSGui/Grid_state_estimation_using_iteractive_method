clc;
clear;
digits(7);
%1)

o1=0;
o2=0;
v1=1;
v2=1;
v3=1;

% define-se o vetor de estados
x=[o1,o2,v1,v2,v3];

syms V1 V2 V3 P31 P32 Q13 Q23 Pg3 Qg3 O2;
%cria-se o vetor de medições
ZOriginal=[1.02,0.9,1.06,-1,0.8,-0.4,-0.2,1.7,0.6,-0.1];
Z=ZOriginal;
Z_simb = [V1 V2 V3 P31 P32 Q13 Q23 Pg3 Qg3 O2];

% W é a matriz dos pesos
W=[1000,0,0,0,0,0,0,0,0,0;0,1000,0,0,0,0,0,0,0,0;0,0,1000,0,0,0,0,0,0,0;0,0,0,100,0,0,0,0,0,0;0,0,0,0,100,0,0,0,0,0;0,0,0,0,0,100,0,0,0,0;0,0,0,0,0,0,100,0,0,0;0,0,0,0,0,0,0,100,0,0;0,0,0,0,0,0,0,0,100,0;0,0,0,0,0,0,0,0,0,10000];

%impedancia nas linhas que será usada para formular o hx 
imp=0.1; 
request=0;
syms V1 V2 V3 O1 O2; 
syms V1 V2 V3 P31 P32 Q13 Q23 Pg3 Qg3 O2;
X=[O1,O2,V1,V2,V3];

P31=V3*V1*sin(0-O1)/imp;
P32=V3*V2*sin(0-O2)/imp;
Q13=(V1*V1-(V1*V3*cos(O1-0)))/imp;
Q23=(V2*V2-(V2*V3*cos(O2-0)))/imp;
Pg3=10*(V3*V1*sin(-O1)+V3*V2*sin(-O2));
Qg3=10*(((V3*V3)-V3*V1*cos(-O1))+((V3*V3)-V3*V2*cos(-O2)));

%Matrix com as expressões de calculo de cada uma das medições
hx=[V1,V2,V3,P31,P32,Q13,Q23,Pg3,Qg3,O2];

%Calculo da jacobiana
Hx=jacobian(hx,X);

%Inicio das varias iterações
i=0;
n_remov=0;
%No primeiro while temos um ciclo que acontece até não ser mais possivel
%retirar medições que estão a aumentar a probabilidade de erro do sistema
while 1
    %1)
    % No segundo ciclo que vai se iniciar a baixo temos o calculo teorico
    % os valores de medição e posteriormente os comparamos com os valores
    % efetivamente medidos. Este ciclo tem fim quando o deltaX que é esta
    % diferença entre os valores do vetor de estado a cada iteração for
    % inferior a 10^-4 (Esse valor é arbitrario e configura uma situação de estabilidade)
    while 1

        o1=x(1);
        o2=x(2);
        v1=x(3);
        v2=x(4);
        v3=x(5);
        % substitui-se os valores dos angulos e tensões nos barramentos na
        % matrix do jacobiano
        H=subs(Hx,{O1,O2,V1,V2,V3},{o1,o2,v1,v2,v3});
        h=subs(hx,{O1,O2,V1,V2,V3},{o1,o2,v1,v2,v3});
        % O r aqui é o residuo (valor das medidas em campo comparadas com as calculadas a partir dos angulos e das tensões)
        
        subs(hx,{O1,O2,V1,V2,V3},{o1,o2,v1,v2,v3});
        r=Z-subs(hx,{O1,O2,V1,V2,V3},{o1,o2,v1,v2,v3});
        r=vpa(double(r'));
        deltaX=(inv(transpose(H)*W*H))*transpose(H)*W*r;
        %Atualiza-se o valor do vetor de estado com os novos valores
        %calculados
        
        x=vpa(x+transpose(deltaX));
        
        %Check contitui a minimização da função objetivo com o objetivo de
        %ter um ideia se o método está a convergir e se os valores estão a
        %ficar mais fiaveis.
       
        check=(double(transpose(r)*W*r));

        %fprintf('\n\nIteração %d\ncheck: %f\nx:%f,%f,%f,%f,%f \n',i,check,x(1),x(2),x(3),x(4),x(5));
        %fprintf('r: ');
        %for j=1:1:length(r)
           %fprintf('%f, ',r(j)); 
        %end
        if abs(max(deltaX))<0.0001
          break;
        end
        i=i+1;
    end
    
    
    
    
    %2)
    %Os graus de liberdade existem para que saibamos se temos mais valores
    %de medições do que icognitas na resolução do nosso sistema.
    degree_of_freedom = length(Z)-length(x);
    grau_confianca= chi2cdf(check,degree_of_freedom);
    
    %3)
    
    % A normalização dos residuos é uma forma de podermos comparar medidas
    % que são de ordem de grandeza diferentes.
    C= inv(W)- H *inv(transpose(H)*(W)*(H))*transpose(H);
    C= double(transpose(diag(C)));
    
    for j = 1:length(r)
        %Obtem-se o vetor dos residuos normalizados
        norm_r(j)=abs(r(j))/sqrt(abs(C(j)));
    end
    
    %for debugging
    if request==2
        norm_r(3)=norm_r(3)+0.03;
    end
    
    
    %4)
    %O elemento cuja eliminação terá maior efeito na diminuição dos
    %resisuos é o elemento que está na posição bad_element_position do
    %vetor dos residuos normalizados
    bad_element_position = find(norm_r == max(norm_r));
    %fprintf("\n\nO elemento ruim é %d, os degraus de liberdade são %d", bad_element_position, degree_of_freedom);
    
    %caso queira a execução do programa ser nenhum output basta substituir
    %o zero a seguir por um numero negativo
    
    if(n_remov == 0)
        %caso queira a execução do programa ser nenhum output basta colovar
        %o
        %prompt = 'Quer ler os valores após remover quantas medidas? ';
        %request = input(prompt);
        
        % De acordo com as contas feitas a retirada da 3 medida é aquelea
        % que fornece resultado com maior grau de confiança, devido a isso
        % decidimos que o processo iria parar antes da remoção da 4 medida
        % mas se desejar basta descomentar o codigo a cima que voltará a
        % mostrar o resultado completo
        request=3; 
    end
    if(degree_of_freedom == 1)
        fprintf("Atenção: Para muitas medidas retiradas as contas podem ter números muito baixos cuja aproximação causa erros\n");
        fprintf("Não foram removidas mais do que %d medidas, as informações a seguir são as referente a remoção número %d\n\n", n_remov+1, n_remov+1);
        break;
    end
    
    if(request == n_remov)
        break;
    end
    
    %A partir desse ponto temos a remoção de dentro de todos os vetores de
    %informação relativa a medição que está mal.
    if degree_of_freedom >1 
        n_remov = n_remov +1;
        Z(bad_element_position) = [];
        Z_simb(bad_element_position) = [];
        r(bad_element_position) = [];
        norm_r(bad_element_position) = [];
        hx(bad_element_position) = [];
        W(bad_element_position,:) = [];
        W(:,bad_element_position) = [];
        H(bad_element_position,:) = [];
        Hx(bad_element_position,:) = [];
        prob_erro = [];
    else
        break;
    end
end

%Insira aqui qualquer medida que deseja obter do nosso programa ou digite
%no terminal após a execução

fprintf("\n\n");
fprintf("As medidas presentes são: ");
Z_simb
fprintf("As cujo valor respectivo em pu é: ");
Z
fprintf("O vetor de estados é: ");
X
fprintf("O vetor com as formulas de calculo de cada medida é: ");
transpose(hx)
fprintf("O valor substituido pelo vetor de estados de hx é: ");
transpose(h)
fprintf("E seus resultados com deltaX menor do que 10^-4 é: ");
x
fprintf("A minimização da função objetivo com deltaX menor do que 10^-4 da: ");
check
fprintf("Grau de confiança: ");
grau_confianca
fprintf("Os residuos de cada medição após a aplicação do método com deltaX menor do que 10^-4 da");
r
fprintf("Matriz da Covariancia");
C
fprintf("A mormalização desses mesmos reisiduos segundo a matriz da covariância");
norm_r
fprintf("A medida a ser eliminada devido o maior erro normalizado de forma a eliminar medidas que estão aferando a estimação de estados é:");
Z_simb(bad_element_position)
fprintf("Número de graus de liberdade restantes:");
degree_of_freedom
fprintf("\n\n");

%Digite aqui qualquer outra medida que deseja obter.


