function [frac_r,fractalD,Frac_type] = modelfitting_fractal(q,I,fractalstart,fractalend)

fractalq = log(q);
fractalI = log(I);
index = find(fractalq >= fractalstart & fractalq <= fractalend);
FITfractalq = fractalq(index);
FITfractalI = fractalI(index);
pdh = polyfit(FITfractalq, FITfractalI, 1);
frac_r = abs(pdh(1));
if frac_r > 0 & frac_r < 3
    fractalD = abs(frac_r);
    Frac_type = 'mass or pore fractal';
elseif  frac_r > 3 & frac_r < 4

    fractalD = 6-abs(frac_r);
    Frac_type = 'surface fractal';
else
    fractalD = NaN;
    Frac_type = NaN;
    print('This SAXS data has no fractal structre');
end

end