function [Pc]=PC(SW,PcS)

Pe=3.5;
swcon=0.2;
lambda=2;

Pc=Pe*PcS*((SW-swcon)/(1-swcon)).^(-1/lambda);
