function Trans=TYWhalf(a,b,P,SW,reservoir,fluid)

if(P(a)>P(b))
    [krw,~]=KR(SW(a));
else
    [krw,~]=KR(SW(b));
end

ky_half=(reservoir.dy+reservoir.dy)/((reservoir.dy/reservoir.Ky(a))+(reservoir.dy/reservoir.Ky(b)));
th_half=2/((1/reservoir.th(a))+(1/reservoir.th(b)));
Trans=6.33e-3*(krw*ky_half*reservoir.dx*th_half)/(fluid.muw*fluid.Bw*reservoir.dy);
