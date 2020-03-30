source('genericFunctions.R')

pairwiseCC <- function(d,i,j){
  dCC <- d[d[[i]] > 0 & d[[j]] > 0 & !is.na(d[[i]]) & !is.na(d[[j]]),]

  cc <- cor(dCC[[i]],dCC[[j]],method='spearman')^2

  ndCC <- length(dCC[[i]])

  d2plot <- data.frame(x=dCC[[i]],y=dCC[[j]])

  gScat <- ggplot(d2plot,aes(x=x,y=y)) +
  geom_point(size=.3,alpha=.3) +
    scale_x_log10() + scale_y_log10() +
    xlab(i) + ylab(j) +
    annotation_logticks(sides='bl',
                        long = unit(.2,'cm'),
                        mid = unit(.1,'cm'),
                        short = unit(.1,'cm'),
                        size = .2) +
    ggtitle(label=bquote("Spearman "*R^2*" = "*.(cc)*" ; N = "*.(ndCC)))

  return(list(cc=cc,
              N=nrow(dCC),
              fig=gScat,
            data=dCC))
}

theme7point()
fDF <- read.table('allData.tab',header=TRUE)
fCol <- which(names(fDF) %in% c('T1_HSstrength_FPM_SSDS','T2_HSstrength_FPM_SSDS','Spo11Oligo_strength','zcwpw1','zcwpw1HS','H3K4me3_b6_kp1_peakStr','h3k4m3Lep','h3k4m3Zyg','h3k4m3MM'))
fCC <- fDF[,fCol]
names(fCC) <- c('SSDS(1)','SSDS(2)','k4m3(12dpp)','spo11','k4m3(Le)', 'k4me(MM)', 'k4me(Zy)','Zcwpw1(Pk)','Zcwpw1(Tot)')

fCCok <- fCC[apply(fCC,1,function(x){if(sum(is.na(x))>0 | sum(x==0) > 0){FALSE}else{TRUE}}),]

ord <- c('SSDS(1)','SSDS(2)','spo11','k4m3(12dpp)','k4m3(Le)','k4me(Zy)','k4me(MM)','Zcwpw1(Pk)', 'Zcwpw1(Tot)')
df4CC <- fCCok[,ord]
ccMat1 <- cor(df4CC,method='spearman')^2

gCC <- ggCorMat(ccMat1,flipIt = TRUE,
                noDiagonal = TRUE,yOnRight = TRUE,keepLeadingZeros = TRUE,decimalPlaces = 2)+
  ggtitle(bquote('CC for intervals in all datasets (N = 1,192; Spearman '*R^2*')'))

ggsave(getIMGname('Zcwpw1_strength_v_other_metrics_at_DSBhotspots','PNG',saveLocation = './'),gCC,width = 3,height = 3)
ggsave(getIMGname('Zcwpw1_strength_v_other_metrics_at_DSBhotspots','PDF',saveLocation = './'),gCC,width = 3,height = 3)

ccMat2 <- ccMat1

for (nI in 1:length(ord)){
  for (nJ in 1:length(ord)){

    pngName <- (getIMGname(paste0(ord[nI],'_v_',ord[nJ],'_strength_Scatter'),'PNG',saveLocation = './'))

    pngName <- gsub("[()]","",gsub("[()]","",pngName))
    pdfName <- gsub('png','pdf',pngName)

    print(pngName)

    lCC <- pairwiseCC(fCC,ord[nI],ord[nJ])
    ccMat2[nI,nJ] <- lCC$cc

    ggsave(pngName,lCC$fig,width = 3,height = 3)
    ggsave(pdfName,lCC$fig,width = 3,height = 3)
  }
}

gCC <- ggCorMat(ccMat2,flipIt = TRUE,
                noDiagonal = TRUE,yOnRight = TRUE,keepLeadingZeros = TRUE,decimalPlaces = 2) +
  ggtitle(bquote('Pairwise CC (N > 2,300; Spearman '*R^2*')'))

ggsave(getIMGname('Zcwpw1_strength_v_other_metrics_at_DSBhotspots_PairwiseCC','PNG',saveLocation = './'),gCC,width = 3,height = 3)
ggsave(getIMGname('Zcwpw1_strength_v_other_metrics_at_DSBhotspots_PairwiseCC','PDF',saveLocation = './'),gCC,width = 3,height = 3)
