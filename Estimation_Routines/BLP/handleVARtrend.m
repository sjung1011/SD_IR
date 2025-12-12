function res=handleVARtrend(Y,Ylag,vecB,modelSpec)
%removes deterministic component from y; expects B=[n(n*nL+1)x1]

nB=numel(vecB); n=modelSpec.modelSize; nL=modelSpec.nVARlags; nH=modelSpec.nHorizons;

T=size(Y,1);

trend=NaN(T+nH+1,n*nL+1); trend(:,1)=1;
trend(1,2:end)=Ylag(1,:);
for j=1:T+nH
    trend(j+1,2:n+1)=trend(j,:)*reshape(vecB,nB/n,n);
    trend(j+1,n+2:end)=trend(j,2:end-n);
end
trend=trend(:,2:n+1);
x=Y-trend(2:T+1,:);

% figure; plotRows=ceil(n/3); pln=1;
% for j=1:n
%     subplot(plotRows,ceil(n/plotRows),pln)
%     
%     plot(Y(:,j)); hold on
%     plot( trend(2:T,j),'--r'); axis tight
%     pln=pln+1;
% end
% set(gcf,'PaperUnits','centimeters','PaperSize',[18 15]) %[x y]
% set(gcf,'PaperPosition',[-1 0 20 15]) %[left bottom width height]
% print(gcf,'-dpdf','VARtrend.pdf'); saveas(gcf,'VARtrend.fig'); 

res.detrended=[zeros(nL,n);x];
res.trend=trend(2:end,:);       %includes horizons
