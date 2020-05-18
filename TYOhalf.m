function Trans=TYOhalf(a,b,P,SW,reservoir,fluid)


if(P(a)>P(b))
    [~,kro]=KR(SW(a));
else
    [~,kro]=KR(SW(b));
end

ky_half=(reservoir.dy+reservoir.dy)/((reservoir.dy/reservoir.Ky(a))+(reservoir.dy/reservoir.Ky(b)));
th_half=2/((1/reservoir.th(a))+(1/reservoir.th(b)));
Trans=6.33e-3*(kro*ky_half*reservoir.dx*th_half)/(fluid.muo(a)*fluid.Bo(a)*reservoir.dy);
