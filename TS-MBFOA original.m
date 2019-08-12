%

clear all
close all
clc

finales=zeros(30,8);  % Matriz para guardar el resultado siguiente en las
% columnas: Mejor, media, Mediana, desv. Estandar, Peor, Tasa de factibilidad,
% Tasa de exito, Rendimiento exitoso de las 30 ejecuciones realiazadas.

corridas=1;  % Ejecuciones independientes [1,30]
problems=100; % Número del problema a resolver (100 es resorte, 101 recipiente)
semillas=zeros(problems,corridas);


for problema=100:problems

    var=numVariables(problema);
    resultado=zeros(30,var+2);      % Matriz que guarda la mejor solucion de cada ejecución.
    sccp=zeros(corridas,1);
    factibles=zeros(corridas,1);
    convergencia15=zeros(1000,2);   % Guardar la convergencia de la ejecucion 15.

    for run=1:corridas

        semilla=sum(100*clock); %% linea del algoritmo
        semillas(problema,run)=semilla; %% linea del algoritmo
        rand('state',semilla); %% linea del algoritmo
        amt0=clock; %% linea del algoritmo
        sb=40;          % Número de población de bacterias [10,500] %% linea del algoritmo
        nc=24;          % Número de ciclos quimiotacticos [1,sb] %% linea del algoritmo
        factorE=1.9;    % Factor de escalamiento [0,2]  %% linea del algoritmo
        repro=1;        % Número de bacterias a reproducir [1, sb/2] %% linea del algoritmo
        repcycle= 100;  % Frecuencia de reproducción [1,gmax/2]  %% linea del algoritmo
        tpasos=0.001;   % [0,1]  %% linea del algoritmo
        eval=15000;     % Número de evaluaciones  [15,000 ->  30,000]   %% linea del algoritmo

        prob=problema;  % Número del problema a resolver (agregarlo a funcion_obj()) %% linea del algoritmo
        contador=0;     % Contador de evaluaciones realizadas por el algoritmo %% linea del algoritmo
        gene=0;     %% linea del algoritmo
        var=numVariables(prob); %% linea del algoritmo
        rango=rangoVariables(prob); %% linea del algoritmo
        valores=eval/(sb*nc); %% linea del algoritmo
        band=1;
        bp=0;
        sperformance=0;
        gr=0;



      if((valores-0.5)>= floor(valores))%% linea del algoritmo
            gmax=floor(valores); %% linea del algoritmo
      else %% linea del algoritmo
          gmax=floor(valores)-1; %% linea del algoritmo
      end %% linea del algoritmo
            vEvals=zeros(gmax,1); %% linea del algoritmo
            nadosExitososG=zeros(gmax,sb);
            primero=0;
            convergencia=zeros(gmax,corridas);
            convergenciaC=zeros(gmax,corridas);
            constraintsvr=zeros(gmax,corridas);
            bacterias=poblacionMecanismoSesgo(sb,prob); % Población aleatoria
            contadorB=0;

     while((contador +(sb*nc)+1)<=eval)   %%codigo del algoritmo %% linea del algoritmo
            gene=gene+1;      %% linea del algoritmo
            bacterias=evaluacions(bacterias,prob); %% linea del algoritmo
            contador=contador+sb; %% linea del algoritmo
            bacteriasSuplente=quimiotaxis(bacterias,prob,nc,factorE);% Ciclo quimiotáctico %% linea del algoritmo
            bacterias(:,1:var)=bacteriasSuplente(:,(1:var));  %% linea del algoritmo
            bacterias(:,(var+var+1):(var+var+2))=bacteriasSuplente(:,(var+var+1:var+var+2)); %% linea del algoritmo

            for h=1:sb
            nadosExitososG(gene,h)=bacteriasSuplente(h,(var+var+3));
            end

            contador=contador+(nc*sb);%% linea del algoritmo

            if(rem(gene,repcycle)==0)  %% linea del algoritmo
                bacterias=reproduccion(bacterias,repro); %Reproducción de bacterias %% linea del algoritmo
            end %% linea del algoritmo

            bacterias=ordenar(bacterias,prob);           %Ordenar bacterias %% linea del algoritmo

            bacterias(sb,:)=eliminacionD(prob);          %Eliminación-D de bacterias %% linea del algoritmo

            contador=contador+1; %% linea del algoritmo

            %%%contadorf=0; non se si se usa esto
       %tamaño de paso estático

bacterias(:,:)=TamPasoEst(prob,tpasos,bacterias); %% linea del algoritmo


       %Tamaño de paso dinamico
%bacterias(:,:)=TamPasoDin(prob,bacterias);  %% linea del algoritmo


           %%%% factibles(run)=contadorf; no se si se usa
%%%%***********metrica de rendimiento exitoso**************
           if(bacterias(1,(var+var+2))==0)
           sperformance=sperformance+succesPerformance(bacterias(1,:),prob); %variable que permite conocer la evaluación donde se cumple la condición de éxito
            if((bp==0) && (sperformance==1))
              sccp(run,1)=contador;
              bp=1;
              contadorMaxB=contador;
            end
           end
%%%%***********metrica de rendimiento exitoso**************

             disp(bacterias(1,var+var+1)); % %% linea del algoritmo visualizacion del resultado alf¿gortimo
            disp(bacterias(1,var+var+2));  % %% linea del algoritmo visualizacion del resultado algoritmo

            bacterias=ordenar(bacterias,prob);

             convergenciaC(gene,1)=bacterias(1,var+var+1);
             convergenciaC(gene,2)=bacterias(1,var+var+2);
             vEvals(gene,1)=contador;

          if(gene>1) %%medir progress ratio para saber el % de mejora de una solucion ...para proxima tesis
             if(convergencia(gene,run)== convergencia(gene-1,run))
                 contadorB=contadorB+1;
             end
          end
           if(run==15)
             convergencia15(gene,1)= bacterias(1,var+var+1);
             convergencia15(gene,2)=bacterias(1,var+var+2);
           end

     end   %% linea del algoritmo end final del ciclo principal while contador de generaciones

     %metricas---se puede enviar a una función----]
      promediofactibles=mean(factibles);
    disp(promediofactibles);
    disp('termino corrida');
    disp(run);
    disp(bacterias(1,(var+var+1):(var+var+2)));
    bacterias(1,var+1);
            for j=1:var
            resultado(run,j)=bacterias(1,j);
            end
            resultado(run,var+1)=bacterias(1,var+var+1);
            resultado(run,var+2)=bacterias(1,var+var+2);
       %metricas---se puede enviar a una función-----]
    end %final de las corridas de ejecuciones independientes

    fi='funcionm3';
        number=problema;
        archivo=strcat(fi,num2str(number));
        savefile = archivo;
        save(savefile, 'resultado'); %para guardar las variables y resultado de la
%función objetivo y obteniendo las estadisticas básicas de los resultados
 tasaExito=succesRate(resultado,run,var,problema);

        estadisticas=zeros(1,2);
        cons=0;
        for d=1:run
            if(resultado(d,var+2)==0)
                cons=cons+1;
                estadisticas(cons,:)=resultado(d,var+1);
            end
        end
     if (cons>0)
        mejor=min(estadisticas);
        peor=max(estadisticas);
        media=mean(estadisticas);
        mediana=median(estadisticas);
        desvEst=std(estadisticas);
    else
        mejor=0;
        peor=0;
        media=0;
        mediana=0;
        desvEst=0;
     end

     %conociendo a succesperformance
     pd=0;
     pscp=0;
     for v=1:run
           if(sccp(v,1)>0)
               pd=pd+1;
               pscp=pscp+sccp(v,1);
           end
     end
     if (pscp>0 && pd>0)
         pscp=pscp/pd;
     else
         pscp=0;
     end
     %conociendo a succesperformance
            finales(problema,1)=mejor(1,1);
            finales(problema,2)=media(1,1);
            finales(problema,3)=mediana(1,1);
            finales(problema,4)=desvEst(1,1);
            finales(problema,5)=peor(1,1);
            if (cons>0)
                finales(problema,6)=(cons*100)/corridas; %tasa de factibilidad
            else
                finales(problema,6)=0;
            end
            if(tasaExito>0)
                finales(problema,7)=(tasaExito*100)/corridas; %tasa de exito
            else
                finales(problema,7)=0;
            end
            if(pscp>0)
                finales(problema,8)=floor(pscp*corridas)/pd; %rendimiento exitoso
            else
                finales(problema,8)=0;
            end
            disp(finales(problema,:));
       disp('problema');
       problema

       [nd,fc]=size(nadosExitososG);

       for v=1:nd
           totale=0;
           for k=1:fc
             totale=totale+ nadosExitososG(v,k);
           end
           Nadostotales(v,1)=totale;
       end
       %successfullswim rate
       swrt=sb*nc*gmax;
       totalnados=0;
       for v=1:nd
          totalnados=totalnados+Nadostotales(v,1);
       end
       successfullswimrate=(totalnados/swrt)*100;
end

        savefile = 'globalm3';
        save(savefile, 'finales'); %para guardar las variables y resultado de la

plot(vEvals(1:gene,1),convergenciaC(1:gene,1),'b--', 'MarkerSize', 3.5);



filename = 'mutmbfoa01.mat';
save(filename)
