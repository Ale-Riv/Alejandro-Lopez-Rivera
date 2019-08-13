clear all
close all
clc

finales=zeros(30,8);
corridas=30;
problems=10;
%tNorm=0; %la normalizacion 1 corresponde a la division sencilla y la 2 corresponde a la sumatoria
semillas=zeros(problems,corridas);
for problema=10:problems

     var=numVariables(problema);
    resultado=zeros(30,var+2);
    sccp=zeros(corridas,1);
    factibles=zeros(corridas,1);
     distG=zeros(1000,corridas);
   convergencia15=zeros(1000,2);

    for run=1:corridas

        semilla=sum(100*clock);
        semillas(problema,run)=semilla;
        rand('state',semilla);
% % % semilla=217017.700000000;
% % % rand('state',semilla);
        amt0=clock;
        sb=20; % número de población de bacterias   probar con40, 50, 100
        nc=10; %número de ciclos quimiotacticos probar con 12 20 25,24
        factorE=1.6; %factor de escalamiento  mayor q 0 menor que 1 parametro sensible 0.8
        repro=5; %numero de bacterias a reproducir
        eval=20000; %numero de evaluaciones  25000, 30000 50000
        eliminacion=1;  %numero de bacterias a eliminar dispersar
        tPasos=0.26;  %parametro sensible
        prob=problema; %numero del problema a resolver (agregarlo a funcion_obj())
        contador=0; %contador de evaluaciones realizadas por el algoritm
        gene=0;
        var=numVariables(prob);
        rang=rangoVariables(prob);
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
                 distancia=zeros(gmax+1,1);
                 bacterias=evaluacions(bacterias,prob);
                %% bacterias=normalizacion(bacterias,tNorm,prob);
%            D=pdist(bacterias(1:sb,1:var),'euclidean');
      %      distancia(1,1)=sum(D);

   while((contador +(sb*nc)+1)<=eval)

            gene=gene+1;

            contador=contador+sb;
            bacterias=ordenar(bacterias,prob);
            %%bacteriasSuplente=quimiotaxis(bacterias,prob,nc,tNorm, factorE);%ciclo quimiotáctico
            bacteriasSuplente=quimiotaxisSinMutacion(bacterias,prob,nc,factorE);%ciclo quimiotáctico

            bacteriasSuplente=evaluacions(bacteriasSuplente,prob);
            bacterias(:,1:var)=bacteriasSuplente(:,(1:var));
            bacterias(:,(var+var+1):(var+var+2))=bacteriasSuplente(:,(var+var+1:var+var+2));
            for h=1:sb
            nadosExitososG(gene,h)=bacteriasSuplente(h,(var+var+3));
            end
            contador=contador+(nc*sb);
            bacterias=ordenar(bacterias,prob);

           % if(rem(gene,100)==0)
                bacterias=reproduccion(bacterias,repro); %reproduccion de bacterias
            %end

    bacterias=ordenar(bacterias,prob);

            bacterias((sb-(eliminacion-1)):sb,:)=eliminacionD(eliminacion,prob); %eliminación-d de bacterias
        contador=contador+eliminacion;

            gene

           %% bacterias=desnormalizacion(bacterias,tNorm,prob);
          %if(gene> 200 && rem(gene,10)==0)
            for u=1:sb
                for k=1:var
                   bacterias(u,var+k)=(rang(k,1) + (rang(k,2) - rang(k,1)) * rand(1))*tPasos; %original
                 %%%  bacterias(u,var+k)=(rang(k,1) + (rang(k,2) - rang(k,1)))*(gene/gmax);
                    %11 de agosto 0.0005
          %          bacterias(u,var+k)=bacterias(u,k+var)-bacterias(u,k+var)*(gene/gmax);
                end
             end
          %end
          bacterias=evaluacions(bacterias,prob);
         %% bacterias=normalizacion(bacterias,tNorm,prob);


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

             convergenciaC(gene,1)=bacterias(1,var+var+1); %fobjetivo
             convergenciaC(gene,2)=bacterias(1,var+var+2);  %suma de violacion de restricciones
             vEvals(gene,1)=contador;

          if(gene>1)
             if(convergenciaC(gene,run)== convergenciaC(gene-1,run))
                 contadorB=contadorB+1;
             end
          end
% %            if(run==15)
% %              convergencia15(gene,1)= bacterias(1,var+var+1);
% %              convergencia15(gene,2)=bacterias(1,var+var+2);
% %           end

%           D=pdist(bacterias(1:sb,1:var),'euclidean');
 %           distancia(gene+1,1)=sum(D);
            %imprimiendo mejores bacterias
          %%  bacterias=desnormalizacion(bacterias,tNorm,prob);

           bacterias=ordenar(bacterias,prob);
            disp(bacterias(1,var+var+1));
            disp(bacterias(1,var+var+2));

          %%  bacterias=normalizacion(bacterias,tNorm,prob);
   end

   if(run==1)
       convergencia15=convergenciaC;
   end
   %% bacterias=desnormalizacion(bacterias,tNorm,prob);
   bacterias=evaluacions(bacterias,prob);
   bacterias=ordenar(bacterias,prob);
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
            distG(1:gmax+1,run)=distancia(:,1);


            %%nados exitosos


    [nd,fc]=size(nadosExitososG);

       for v=1:nd
           totale=0;
           for k=1:fc
             totale=totale+ nadosExitososG(v,k);
           end
           Nadostotales(v,1)=totale;
       end
       if(run==15)
       NadosTotalesCorrida15=Nadostotales;
       end


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
   %%%de aqui quite nados exitosos

       %successfullswim rate
       swrt=sb*nc*gmax;
       totalnados=0;
       for v=1:nd
          totalnados=totalnados+Nadostotales(v,1);
       end
       if(totalnados>0)
       finales(problema,9)=totalnados; % nados totoales
       else
         finales(problema,9)=0;
       end
       successfullswimrate=(totalnados/swrt)*100; %nados exitosos

        if(totalnados>0)
       finales(problema,10)=successfullswimrate;
        else
         finales(problema,10)=0;
        end

    savefile = 'globalm3n';
        save(savefile, 'finales');

% % if problems==10 && tNorm==0
% %  xlswrite('resultadosNormalizados30.xlsx', resultado, 'Resorte');
% %   xlswrite('resultadosNormalizados30.xlsx', finales(10,:), 'Resorte','H4:Q4');
% % filename = 'mbfoaOriginalResorte30n.mat';
% % end
% %
% %
% % if problems==11 && tNorm==0
% %  xlswrite('resultadosNormalizados30.xlsx', resultado, 'Recipiente');
% %   xlswrite('resultadosNormalizados30.xlsx', finales(11,:), 'Recipiente','H4:Q4');
% % filename = 'mbfoaOriginalRecipiente30n.mat';
% % end
% %
% %
% %
% % if problems==12 && tNorm==0
% %  xlswrite('resultadosNormalizados30.xlsx', resultado, 'Viga');
% %   xlswrite('resultadosNormalizados30.xlsx', finales(12,:), 'Viga','H4:Q4');
% % filename = 'mbfoaOriginalViga30n.mat';
% % end


end

%gdf
        savefile = 'globalm3';
        save(savefile, 'finales'); %para guardar las variables y resultado de la
%función objetivo y
 %disp(bacterias(1,:)); %al final del ejecución se presenta la mejor bacteria
%plot(convergencia(1:gene,1));
% % % for(jh=1:gmax)
% % %  vv(jh,1)=jh;
% % % end
% % % plot(vv(:,1),NadosTotalesCorrida15(:,1),'ro', 'MarkerSize', 3.5);

%plot(vv(:,1),convergencia(:,1),'ro', 'MarkerSize', 3.5);

%%%%%plot(vEvals(:,1),convergenciaC(:,1),'ro', 'MarkerSize', 3.5);
 % para cuando haya solo una ejecucion usar laura
 fig=figure;
plot(vEvals(1:gene,1),convergencia15(1:gene,1),'ro', 'MarkerSize', 3.5); %para cuando
%haya 30 ejecuciones usar laura

%savefile = 'test.mat';
%save(savefile, 'resultado'); %para guardar las variables y resultado de la función objetivo y
%violación de restricciones el resultado en un archivo

% plot(vEvals(1:gene,1),convergenciaC(1:gene,1),'b--', 'MarkerSize', 3.5);
print(fig,'convergencia15','-dpng')

    savefig('convergencia15.fig')


%dfdfdf

filename = 'mbfoa160519.mat';
save(filename)
