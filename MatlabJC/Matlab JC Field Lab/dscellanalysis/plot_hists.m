function[] = plot_hists(mag, dsindex, magmax, magave)

subplot(2,4,1);
hist(mag{1,1},100);
title('Param mag speed 32');
subplot(2,4,2)
hist(mag{2,1},100);
title('Param mag speed 256');
subplot(2,4,3)
plot(mag{1,1}, mag{2,1},'o');
title('Speed 32 vs speed 256 param');
subplot(2,4,4)
plot(log(mag{1,1}), log(mag{2,1}),'o');
title('Speed 32 vs speed 256 log param');

subplot(2,4,5)
hist(magave,100)
title('Ave param');
subplot(2,4,6)
hist(magmax,100)
title('Max Param');

subplot(2,4,7)
hist(dsindex{1,1},100);
title('DS Index speed 32');
subplot(2,4,8)
hist(dsindex{2,1},100);
title('DS Index speed 256');
end

