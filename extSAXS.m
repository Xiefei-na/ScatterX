function [qfull,Ifull,Rg] = extSAXS(q,I,Porodlimit,Guinierlimit,Exlimit)
porodq = q.^2;
porodI = log(q.^4.*I);
index = find(porodq >= Porodlimit & porodq <= porodq(end));
FITporodq = porodq(index);
FITporodI = porodI(index);
pr = polyfit(FITporodq, FITporodI, 1);
caliI = exp( real(porodI - pr(1)*porodq))./(q.^4);
qkuo = (max(q) + (q(3) - q(2)):q(3) - q(2):Exlimit)';
IPorodkuo = exp(pr(2))./(qkuo.^4);
q_right = [q; qkuo];
I_right = [caliI; IPorodkuo];
index2 = find(q_right >= q_right(1) & q_right <= Guinierlimit);
qc = q_right(index2);
Ic = I_right(index2);
p = polyfit(qc.^2, log(Ic), 1);

Rg = sqrt(-3*p(1));
qguinier = 0:0.01:q_right(1);
Iguinier = exp(p(1).*qguinier.^2+p(2));
qfull = [qguinier';q_right];
Ifull = [Iguinier';I_right];
end