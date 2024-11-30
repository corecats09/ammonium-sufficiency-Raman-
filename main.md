
```matlab
plot (x,mean(log_wt,2))
hold on
plot (x,mean(log_glnb,2))
plot (x,mean(log_glnk,2))
xlabel('wavenumber/cm^{-1}')
ylabel('intensity/counts')
legend('wt','\Deltaglnb','\Deltaglnk')
legend('wt(84)','\Deltaglnb(83)','\Deltaglnk(78)')
% mean spectrum of the samples(number),
title('in the condition of ammonium sufficiency, exponential growth')
```
![meanSpectrum_log_wt_glnb_glnk](figures/mean1.jpg)

```matlab
plot (x,mean(ro_wt,2))
hold on
plot (x,mean(ro_glnb,2))
plot (x,mean(ro_glnk,2))
xlabel('wavenumber/cm^{-1}')
ylabel('intensity/counts')
legend('wt','\Deltaglnb','\Deltaglnk')
legend('wt(76)','\Deltaglnb(88)','\Deltaglnk(88)')
% mean spectrum of the samples(number),
title('in the condition of ammonium running out')
```
![meanSpectrum_log_wt_glnb_glnk](figures/mean2.jpg)
```matlab
plot (x,mean(sta_wt,2))
hold on
plot (x,mean(sta_glnb,2))
plot (x,mean(sta_glnk,2))
xlabel('wavenumber/cm^{-1}')
ylabel('intensity/counts')
legend('wt','\Deltaglnb','\Deltaglnk')
legend('wt(90)','\Deltaglnb(71)','\Deltaglnk(85)')
% mean spectrum of the samples(number),
title('in the condition of ammonium starvation')

```
![meanSpectrum_log_wt_glnb_glnk](figures/mean3.jpg)

```matlab
plot (x,mean(upsft_wt,2))
hold on
plot (x,mean(upsft_glnb,2))
plot (x,mean(upsft_glnk,2))
xlabel('wavenumber/cm^{-1}')
ylabel('intensity/counts')
legend('wt','\Deltaglnb','\Deltaglnk')
legend('wt(85)','\Deltaglnb(80)','\Deltaglnk(84)')
% mean spectrum of the samples(number),
title('in the condition of ammonium up-shift')
```
![meanSpectrum_log_wt_glnb_glnk](figures/mean4.jpg)

```matlab
datasetall = [log_wt,log_glnb,log_glnk,ro_wt,ro_glnb,ro_glnk,sta_wt,sta_glnb,sta_glnk,upsft_wt,upsft_glnb,upsft_glnk];
quantiles = [];
numbers = [];
numbers(1) = width(log_wt);
numbers(2) = width(log_glnb);
numbers(3) = width(log_glnk);
numbers(4) = width(ro_wt);
numbers(5) = width(ro_glnb);
numbers(6) = width(ro_glnk);
numbers(7) = width(sta_wt);
numbers(8) = width(sta_glnb);
numbers(9) = width(sta_glnk);
numbers(10) = width(upsft_wt);
numbers(11) = width(upsft_glnb);
numbers(12) = width(upsft_glnk);
quantiles(1) =0;

for n = 1:12
quantiles(n+1) = sum(numbers(1:n));
end
[coeff,score,latent] = pca(datasetall');

for n = 1:3
scatter3 (score(quantiles(n)+1:quantiles(n+1),1),score(quantiles(n)+1:quantiles(n+1),2),score(quantiles(n)+1:quantiles(n+1),3),12,'o','filled')
hold on
end
```

```matlab
for n = 4:6
scatter3 (score(quantiles(n)+1:quantiles(n+1),1),score(quantiles(n)+1:quantiles(n+1),2),score(quantiles(n)+1:quantiles(n+1),3),12,'^')
hold on
end

for n = 7:9
scatter3 (score(quantiles(n)+1:quantiles(n+1),1),score(quantiles(n)+1:quantiles(n+1),2),score(quantiles(n)+1:quantiles(n+1),3),12,'x')
hold on
end

for n = 10:12
scatter3 (score(quantiles(n)+1:quantiles(n+1),1),score(quantiles(n)+1:quantiles(n+1),2),score(quantiles(n)+1:quantiles(n+1),3),12,'square','filled')
hold on
end





