function fitness_val=fitness(P1,Position,Eo,Sensing_range)
global distance x1 x2 x3
alpha1=0.5;
alpha2=0.2;
alpha3=1-alpha1-alpha2;
for ii=1:size(Position,2)
    distance(ii)=sqrt(((P1(1)-Position(1,ii)).^2)+((P1(2)-Position(2,ii)).^2));
    if distance(ii)<Sensing_range
        min_id(ii)=ii;
    end
end
min_id(min_id==0)=[];
Cn=length(min_id);

x1=sum(distance)/Cn;
x2=((Eo*Cn)/Cn)/Eo;
x3=1/Cn;

fitness_val=(alpha1*x1)+(alpha2*x2)+(alpha3*x3);