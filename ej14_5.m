ro = 0.294;
km = 10000;
xmin = 0;
xmax = 12;
ymin = 1000;

yex = (km * ymin * exp(ro * xmax)) / (km + ymin * (exp(ro * xmax) -1));

for NI = 1:100

    h = (xmax - xmin) / (NI + 1);
    x = linspace(xmin, xmax, NI + 2)';
    ysime = [ymin];

    for i = 1:length(x)-1
        yip1 = (ysime(end) * (ro * h + 1)) - (ro * h * (ysime(end)^2) / km);
        ysime = [ysime; yip1];
    endfor

    err = 100 * (abs(ysime(end) - yex) / yex);

    if err < 1
        err
        NI
        h
        break;
    end

end
