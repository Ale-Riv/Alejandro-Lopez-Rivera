function bact=reproduccion(bacteria, Sr)

[f,c]=size(bacteria);

   for y=(f-Sr)+1:f
       for h=1:c
           bacteria(y,h)=bacteria((f-y)+1,h);
       end
   end
bact=bacteria;
