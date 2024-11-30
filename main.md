
```matlab
plot (x,mean(log_wt,2))
hold on
plot (x,mean(log_glnb,2))
plot (x,mean(log_glnk,2))
xlabel('wavenumber/cm^{-1}')
ylabel('intensity/counts')
legend('wt','\Deltaglnb','\Deltaglnk')
legend('wt(84)','\Deltaglnb(83)','\Deltaglnk(78)')
```
![meanSpectrum_log_wt_glnb_glnk](figures/mean1.jpg)
