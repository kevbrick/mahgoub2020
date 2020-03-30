source('genericFunctions.R')

theme7point()

zD <- read.table('allData.tab',header=TRUE)

combos <- c('Spo11Oligo_strength',
            'T1_HSstrength_FPM_SSDS',
            'O1_HSstrength_FPM_SSDS',
            'H3K4me3_b6_kp1_peakStr',
            'h3k4m3Zyg',
            'Hop2ko_HSstrength_FPM_SSDS',
            'prdm9',
            'h3k36m3',
            'h3k36m3SPO')

names  <- c('Spo11 Oligos',
            'DMC1 SSDS (testis)',
            'DMC1 SSDS (ovary)',
            'Hotspot H3K4me3 (12dpp)',
            'Hotspot H3K4me3 (Zygonema)',
            'SSDS Hop2ko',
            'prdm9',
            'h3k36m3',
            'h3k36m3SPO')

combosS <- c('zcwpw1HS',
            'T1_HSstrength_FPM_SSDS',
            'O1_HSstrength_FPM_SSDS',
            'H3K4me3_b6_kp1_peakStr',
            'h3k4m3Zyg',
            'Hop2ko_HSstrength_FPM_SSDS',
            'prdm9',
            'h3k36m3',
            'h3k36m3SPO')

namesS  <- c('Zcwpw1',
            'DMC1 SSDS (testis)',
            'DMC1 SSDS (ovary)',
            'Hotspot H3K4me3 (12dpp)',
            'Hotspot H3K4me3 (Zygonema)',
            'SSDS Hop2ko',
            'prdm9',
            'h3k36m3',
            'h3k36m3SPO')

dList <- list()
sList <- list()

zPlot <- zD[zD$zcwpw1>0 & !zD$DefaultHS_Oltotal,]
for (c in 1:length(combos)){
  dList[[combos[c]]] <- ggSmoothScatter(zPlot,
                                        col1 = 'zcwpw1HS',
                                        col2 = combos[[c]],
                                        useLog = TRUE,
                                        nRes = 35,
                                        xlab = 'Zcwpw1 signal',ylab=names[c],
                                        midCol='black',omitZeros = TRUE,forScreen = FALSE)

  sList[[combosS[c]]] <- ggSmoothScatter(zPlot,
                                        col1 = 'Spo11Oligo_strength',
                                        col2 = combosS[[c]],
                                        useLog = TRUE,
                                        nRes = 35,
                                        xlab = 'Spo11 signal',ylab=namesS[c],
                                        midCol='black',omitZeros = TRUE,forScreen = FALSE)
}

gZ <- ggarrange(dList$Spo11Oligo_strength$fig,
          dList$T1_HSstrength_FPM_SSDS$fig,
          dList$O1_HSstrength_FPM_SSDS$fig,
          dList$H3K4me3_b6_kp1_peakStr$fig,
          dList$h3k4m3Zyg$fig,
          dList$Hop2ko_HSstrength_FPM_SSDS$fig,
          dList$prdm9$fig,
          dList$h3k36m3$fig,
          dList$h3k36m3SPO$fig,
          ncol=3,nrow=3)

ggsave('scatters_Zcw_v_recombination.pdf',height=5,width=5)
ggsave('scatters_Zcw_v_recombination.png',height=5,width=5,dpi=400,units='in')

gS <- ggarrange(sList$zcwpw1HS$fig,
          sList$T1_HSstrength_FPM_SSDS$fig,
          sList$O1_HSstrength_FPM_SSDS$fig,
          sList$H3K4me3_b6_kp1_peakStr$fig,
          sList$h3k4m3Zyg$fig,
          sList$Hop2ko_HSstrength_FPM_SSDS$fig,
          sList$prdm9$fig,
          sList$h3k36m3$fig,
          sList$h3k36m3SPO$fig,
          ncol=3,nrow=3)

ggsave('scatters_Spo11_v_recombination.pdf',height=5,width=5)
ggsave('scatters_Spo11_v_recombination.png',height=5,width=5,dpi=400,units='in')
