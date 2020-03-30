source('genericFunctions.R')

theme7point()

c <- fread('allSPoTVals.tab')
cM <- reshape2:::melt.data.frame(c,id.vars = c("bam","interval"),measure.vars=c("pcInside"))

zV <- cM[cM$bam == 'GFP_Control_B6',c('interval','value')]
names(zV) <- c('interval','GFP')
cMPC <- plyr:::join(cM,zV,by='interval')

gA <- ggplot(cM,aes(x=bam,y=interval,fill=value,label=paste0(round(value,1),'%'))) + 
  geom_tile() + 
  scale_fill_gradient2("SPoT (%)",low='white',high='red') + 
  geom_text(size=7*5/14) + 
  geom_vline(xintercept=2.5) + 
  xlab('Experimental dataset') + 
  ylab('Genomic intervals') + 
  theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1),
        legend.position='top',
        legend.direction = 'horizontal',
        legend.key.width = unit(1,'cm'))

gB <- ggplot(cMPC,aes(x=bam,y=interval,fill=log(value/GFP),label=round(value/GFP,1))) + 
  geom_tile() + 
  scale_fill_gradient2(bquote(log*"(signal/"*GFP[B6]*")"),
                       low='dodgerblue3',mid='white',high='red',
                       midpoint=log(1)) + 
  geom_text(size=7*5/14) + 
  geom_vline(xintercept=2.5) + 
  xlab('Experimental dataset') + 
  ylab('Genomic intervals') + 
  theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1),
        legend.position='top',
        legend.direction = 'horizontal',
        legend.key.width = unit(1,'cm'))

gX2 <- ggarrange(gA,gB,
                 ncol=1,nrow=2,
                 labels=c('a','b'),
                 font.label = list(size=8,face='bold'))

svg('SPoT_Tables.svg',width=4,height=7)
print(gX2)
dev.off()

pdf('SPoT_Tables.pdf',width=4,height=7)
print(gX2)
dev.off()

png('SPoT_Tables.png',width=4,height=7,res=400,units='in')
print(gX2)
dev.off()