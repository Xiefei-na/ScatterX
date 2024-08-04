function [Fq,FI] = cutuseless(q,inten,q1,q2)

index = find(q >= q1 & q <= q2);
Fq = q(index);
FI = inten(index);

end
