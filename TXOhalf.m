function Trans=TXOhalf(a,b,P,SW,reservoir,fluid)

if(P(a)>P(b))
    [~,kro]=KR(SW(a));
else
    [~,kro]=KR(SW(b));
end
    
kx_half=(reservoir.dx+reservoir.dx)/((reservoir.dx/reservoir.Kx(a))+(reservoir.dx/reservoir.Kx(b)));
th_half=2/((1/reservoir.th(a))+(1/reservoir.th(b)));
Trans=6.33e-3*(kro*kx_half*reservoir.dy*th_half)/(fluid.muo(a)*fluid.Bo(a)*reservoir.dx);
