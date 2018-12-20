function Error = PairCorrAnalyticalSolutionfit(a,Ve1e2,Vi1i2,Ve1i2,Vi1e2,Ve1,Vi1,Ve1i1,Ve2,Vi2,Ve2i2,Pss) 

numerator = Ve1e2 + a^2*Vi1i2 - a*Ve1i2 - a*Vi1e2 ;

denominator = sqrt((Ve1 + a^2*Vi1 - 2*a*Ve1i1)*(Ve2 + a^2*Vi2 - 2*a*Ve2i2)) ;

Error = numerator/denominator - Pss ;

end