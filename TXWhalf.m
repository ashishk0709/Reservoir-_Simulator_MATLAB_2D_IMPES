function Trans=TXWhalf(a,b,P,SW,reservoir,fluid)

if(P(a)>P(b))
    [krw,~]=KR(SW(a));
else
    [krw,~]=KR(SW(b));
end
    
kx_half=(reservoir.dx+reservoir.dx)/((reservoir.dx/reservoir.Kx(a))+(reservoir.dx/reservoir.Kx(b)));
th_half=2/((1/reservoir.th(a))+(1/reservoir.th(b)));
Trans=6.33e-3*(krw*kx_half*reservoir.dy*th_half)/(fluid.muw*fluid.Bw*reservoir.dx);
