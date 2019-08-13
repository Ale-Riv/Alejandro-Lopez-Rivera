clear all
close all
clc

     semilla=sum(100*clock);
     rand('state',semilla);
     amt0=clock;

      %*************PROBLEMA A RESOLVER****************

        Var=3;          % Número de variables del problema a resolver (Resorte de tensión/compresión)

        Rang=[0.05 2;     %X1 -- X=(0.05,2)
              0.25 1.3;   %X2 -- Y=(0.25,1.23)
              2 15];      %X3 -- Z=(2,15)  ************** Rango de variables del problema

        %*************PROBLEMA A RESOLVER***************

        %*************PARÁMETROS DE ENTRADA***************

        Sb=20;          % Número de población de bacterias [10,500]
        Nc=10;          % Número de ciclos quimiotáxicos [1,Sb]
        B=1.6;          % Factor de escalamiento **BETA** [0,2]
        Sr=5;           % Número de bacterias a reproducir [1, Sb/2]
        R=0.026;        % Tamaño de paso estático[0,1]
        eval=20000;     % Número de evaluaciones  [15,000 ->  30,000]

        %*************PARÁMETROS DE ENTRADA***************

        contador=0;     % Contador de evaluaciones realizadas por el algoritmo
        gene=0;
        valores=eval/(Sb*Nc);

        corridas=1;
        resultado=zeros(30,Var+2);  % Matriz que guarda la mejor solución de cada ejecución.
        factibles=zeros(corridas,1);
        run=1:corridas;

      if((valores-0.5)>= floor(valores))
          gmax=floor(valores);
         else
          gmax=floor(valores)-1;
      end
      if(gmax<50)
       RepCycle=2;
      else
          RepCycle=20;
      end
            vEvals=zeros(gmax,1);
            nadosExitososG=zeros(gmax,Sb);

            bacterias=poblacionMecanismoSesgo(Sb,Rang,Var);                % Población aleatoria

            contadorB=0;

     while((contador +(Sb*Nc)+1)<=eval)   %Inicia procesos posteriores a la generación de la población
            gene=gene+1;
            disp(gene);

            bacterias=evaluacions(bacterias,Var);                          %Evaluaciones

            contador=contador+Sb;

            bacteriasSuplente=quimiotaxis(bacterias,Nc,B,Var,Rang);        % Ciclo quimiotáxico

            bacterias(:,1:Var)=bacteriasSuplente(:,(1:Var));
            bacterias(:,(Var+Var+1):(Var+Var+2))=bacteriasSuplente(:,(Var+Var+1:Var+Var+2));

            for h=1:Sb
              nadosExitososG(gene,h)=bacteriasSuplente(h,(Var+Var+3));
            end

            contador=contador+(Nc*Sb);

            if(rem(gene,RepCycle)==0)
              bacterias=reproduccion(bacterias,Sr);                        % Reproducción de bacterias
            end

            bacterias=ordenar(bacterias,Var);                              % Ordenar bacterias

            bacterias(Sb,:)=eliminacionD(Var,Rang);                        % Eliminación dispersión

            contador=contador+1;

            bacterias(:,:)=TamPasoEst(Rang,R,bacterias,Var);               %Tamaño de paso estático

            %bacterias(:,:)=TamPasoDin(rango,bacterias,Var);               %Tamaño de paso dinámico


            disp(bacterias(1,Var+Var+1)); % Visualización de resultado del algoritmo
            disp(bacterias(1,Var+Var+2)); % Visualización de resultado del algoritmo

           vEvals(gene,1)=contador;

     end   % Final del ciclo principal (while contador de generaciones)


    promediofactibles=mean(factibles);
    disp(promediofactibles);
    disp('termino corrida');
    disp(run);
    disp(bacterias(1,(Var+Var+1):(Var+Var+2)));
    bacterias(1,Var+1);
      for j=1:Var
            resultado(run,j)=bacterias(1,j);
      end
            resultado(run,Var+1)=bacterias(1,Var+Var+1);
            resultado(run,Var+2)=bacterias(1,Var+Var+2);

filename = 'mutmbfoa01.mat';
save(filename)
