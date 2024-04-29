a<-read.csv('caribbean_metrics.csv',header=T)

ok<-which(is.finite(rowSums(a)))
a<-a[ok,]

a<-data.frame(a)

pdf('new_weighted_metrics_ok.pdf',height=8,width=8)
par(mfrow=c(2,2))

r2<-round(cor.origin(a$II,a$I_c)**2,2)
plot(a$II,a$I_c,las=1,cex.axis=1.2,
     cex.lab=1.2,pch=16,xlab='mean prevalence',ylab='weighted mean prevalence',main=r2)

abline(lm(a$II~a$I_c),col='red')
abline(lm(a$II~0+a$I_c))
abline(0,1,lty=2)


r2<-round(cor.origin(a$I,a$wI)**2,2)
plot(a$I,a$wI,las=1,cex.axis=1.2,
     cex.lab=1.2,pch=16,xlab='cumulative prevalence',ylab='weighted cumulative prevalence',main=r2)

abline(lm(a$I~a$wI),col='red')
abline(lm(a$I~0+a$wI))
abline(0,1,lty=2)


r2<-round(cor.origin(a$IV,a$III_c)**2,2)
plot(a$IV,a$III_c,las=1,cex.axis=1.2,
     cex.lab=1.2,pch=16,xlab='mean mortality',ylab='weighted mean mortality',main=r2)

abline(lm(a$IV~a$III_c),col='red')
abline(lm(a$IV~0+a$III_c))
abline(0,1,lty=2)

r2<-round(cor.origin(a$III,a$wIII)**2,2)
plot(a$III,a$wIII,las=1,cex.axis=1.2,
     cex.lab=1.2,pch=16,xlab='cumulative mortality',ylab='weighted cumulative mortality',main=r2)

abline(lm(a$III~a$wIII),col='red')
abline(lm(a$III~0+a$wIII))
abline(0,1,lty=2)


dev.off()


colnames(a)<-c('cumulative prevalence',
               'weighted mean prevalence',
               'mean prevalence',
               'II_c',
               'cumulative mortality',
               'weighted mean mortality',
               'mean mortality',
               'IV_c',
               'weighted cumulative prevalence',
               'weighted cumulative mortality',
               'coral_genera_n',
               'total_coral_cover',
               'coral_cover_sd')


sink('coefficients_w_vs_u.txt')
print(lm(a$`weighted mean prevalence`~a$`mean prevalence`))
print(lm(a$`weighted cumulative prevalence`~a$`cumulative prevalence`))
print(lm(a$`weighted mean mortality`~a$`mean mortality`))
print(lm(a$`weighted cumulative mortality`~a$`cumulative mortality`))

print('forced to the origin')
print(lm(a$`weighted mean prevalence`~0+a$`mean prevalence`))
print(lm(a$`weighted cumulative prevalence`~0+a$`cumulative prevalence`))
print(lm(a$`weighted mean mortality`~0+a$`mean mortality`))
print(lm(a$`weighted cumulative mortality`~0+a$`cumulative mortality`))

sink()

write.table(a,'data_weighted_vs_unweighted.csv',quote=F,sep=',',
            row.names=F,col.names=T)


