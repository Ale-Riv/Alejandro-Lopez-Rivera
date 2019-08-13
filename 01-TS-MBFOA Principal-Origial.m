clear all
close all
clc

finales=zeros(30,8);  %Matriz para guardar el resultados de las ejecuciones independientes con las siguientes columnas
% Mejor, media, mediana, desv. Estandar, peor, tasa de factibilidad, tasa
% de exito, rendimiento exitoso.
corridas=30;  %ejecuciones independientes [1,30]
problems=100; %numero del problema a resolver
semillas=zeros(problems,corridas);
for problema=100:problems

     var=numVariables(problema);
    resultado=zeros(30,var+2);
    sccp=zeros(corridas,1); %
    factibles=zeros(corridas,1);  %
   convergencia15=zeros(1000,2);  %guardar la convergencia de la ejecucion 15

    for run=1:corridas

        semilla=sum(100*clock);
        semillas(problema,run)=semilla;
rand('state',semilla);
        amt0=clock;
        sb=40; % número de población de bacterias [10,500]
        nc=24; %número de ciclos quimiotacticos [1,sb]
        factorE=1.9; %factor de escalamiento [0,2]
        repro=1; %numero de bacterias a reproducir [1, sb/2]
        repcycle= 100; %frecuencia de reproducción [1,gmax/2]
        tpasos=0.001; %[0,1]
        eval=30000; %numero de evaluaciones  [15,000 ->  30,000]

        prob=problema; %numero del problema a resolver (agregarlo a funcion_obj())
        contador=0; %contador de evaluaciones realizadas por el algoritm
        gene=0;
        var=numVariables(prob);
        rango=rangoVariables(prob);
        valores=eval/(sb*nc);
        band=1;
        bp=0;
        sperformance=0;
        gr=0;
      if((valores-0.5)>= floor(valores))
            gmax=floor(valores);
      else
          gmax=floor(valores)-1;
      end
        vEvals=zeros(gmax,1);
        nadosExitososG=zeros(gmax,sb);
        primero=0;
         convergencia=zeros(gmax,corridas);
         convergenciaC=zeros(gmax,corridas);
          constraintsvr=zeros(gmax,corridas);
         bacterias=poblacion(sb,prob); %poblacion aleatoria
                contadorB=0;

   while((contador +(sb*nc)+1)<=eval)
            gene=gene+1;
            bacterias=evaluacions(bacterias,prob);
            contador=contador+sb;
            bacteriasSuplente=quimiotaxis(bacterias,prob,nc);%ciclo quimiotáctico
            bacterias(:,1:var)=bacteriasSuplente(:,(1:var));
            bacterias(:,(var+var+1):(var+var+2))=bacteriasSuplente(:,(var+var+1:var+var+2));

            for h=1:sb
            nadosExitososG(gene,h)=bacteriasSuplente(h,(var+var+3));
            end

            contador=contador+(nc*sb);
            if(rem(gene,repcycle)==0)
                bacterias=reproduccion(bacterias,repro); %reproduccion de bacterias
            end

            bacterias=ordenar(bacterias,prob);

            bacterias(sb,:)=eliminacionD(prob); %eliminación-d de bacterias
            contador=contador+1;
            gene
            %******tamaño de paso estático ******
          for u=1:sb
                for k=1:var
                    vector(u,var+k)=(rango(k,1) + (rango(k,2) - rango(k,1)))*tpasos;
                end
          end
            %******tamaño de paso estático ******
           % for u=1:sb
           %     for k=1:var
           %         vector(u,var+k)= bacterias()*(gene/gmax);
           %     end
          %end

            %****tamaño de paso dinámico

            contadorf=0;
           for d=1:sb
              if(bacterias(d,var+var+2)==0)
                  contadorf=contadorf+1;
              end
           end
             factibles(run)=contadorf;

           if(bacterias(1,(var+var+2))==0)
           sperformance=sperformance+succesPerformance(bacterias(1,:),prob); %variable que permite conocer la evaluación donde se cumple la condición de éxito
            if((bp==0) && (sperformance==1))
              sccp(run,1)=contador;
              bp=1;
              contadorMaxB=contador;
            end
           end

             disp(bacterias(1,var+var+1));
            disp(bacterias(1,var+var+2));
             bacterias=ordenar(bacterias,prob);

             convergenciaC(gene,1)=bacterias(1,var+var+1);
             convergenciaC(gene,2)=bacterias(1,var+var+2);
             vEvals(gene,1)=contador;

          if(gene>1)
             if(convergencia(gene,run)== convergencia(gene-1,run))
                 contadorB=contadorB+1;
             end
          end
           if(run==15)
             convergencia15(gene,1)= bacterias(1,var+var+1);
             convergencia15(gene,2)=bacterias(1,var+var+2);
           end

      end
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

end %final de las corridas de ejecuciones independientes

    fi='funcionm3';
        number=problema;
        archivo=strcat(fi,num2str(number));
        savefile = archivo;
        save(savefile, 'resultado'); %para guardar las variables y resultado de la
%función objetivo y obteniendo las estadisticas básicas de los resultados
 tasaExito=succesRate(resultado,run,var,problema);

   % if ((problema == 1) || (problema ==2) || (problema == 3) || (problema == 4) || (problema ==5) || (problema == 6) || (problema == 7) || (problema ==8) || (problema == 9) || (problema == 10) || (problema ==11) || (problema == 12) || (problema == 13) || (problema ==15) || (problema == 16) || (problema == 17) || (problema ==18) || (problema == 19) || (problema == 20) || (problema ==23) || (problema == 24))
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
%función objetivo y
 %disp(bacterias(1,:)); %al final del ejecución se presenta la mejor bacteria
%plot(convergencia(1:gene,1));
%for(jh=1:gmax)
 %vv(jh,1)=jh;
%end
%plot(vv(:,1),Nadostotales(:,1),'ro', 'MarkerSize', 3.5);

%plot(vv(:,1),convergencia(:,1),'ro', 'MarkerSize', 3.5);

plot(vEvals(:,1),convergenciaC(:,1),'ro', 'MarkerSize', 3.5);

%savefile = 'test.mat';
%save(savefile, 'resultado'); %para guardar las variables y resultado de la función objetivo y
%violación de restricciones el resultado en un archivo


filename = 'mutmbfoa01.mat';
save(filename)
