################# ################# ################# ################# 
################# Other Similarity calculation method ##############

#########Pearson  ##################################
dataraw = read.csv(".../data/dnapro_f/phage_dna_pro.csv",header = F)
dataraw = t(dataraw)
data1 = cor(dataraw)
write.csv(data1,file = ".../data/dnapro_f/phage_dnapro_corsim.csv",row.names = F)

#########spearman  ##################################
dataraw = read.csv(".../data/dnapro_f/phage_dna_pro.csv",header = F)
dataraw = t(dataraw)
data1 = cor(dataraw, method = "spearman" )
write.csv(data1,file = ".../data/dnapro_f/phage_dnapro_spearmansim.csv",row.names = F)

#########kendall  ##################################
dataraw = read.csv(".../data/dnapro_f/phage_dna_pro.csv",header = F)
dataraw = t(dataraw)
data1 = cor(dataraw, method = "kendall" )
write.csv(data1,file = ".../data/dnapro_f/phage_dnapro_kendallsim.csv",row.names = F)

######### manhattan ##################################
dataraw = read.csv(".../data/dnapro_f/phage_dna_pro.csv",header = F)
dataraw = t(dataraw)
data1 = as.data.frame(dataraw)
data1 = apply(data1, 1, as.numeric)
d = dist(data1, method = "manhattan",diag = T,upper = T)
d = as.matrix(d)
esim = 1/(1+d)     
write.csv(esim,file = ".../data/dnapro_f/phage_dnapro_manhattansim.csv",row.names = F)




################# ################# ################# ################# 
#################    KNN graph   ##################
data = read.csv(".../data/dnapro_f/phage_dnapro_sim.csv",header = F)
data1 = as.data.frame(data)
data1 = apply(data1, 1, as.numeric)
K = 6  ##the number of adaptive neighbors
for (i in 1:nrow(data1)){
  datai = as.numeric(data1[i,])
  datai[which(datai < sort(unique(datai),decreasing = T)[K])]=0
  data1[i,] = datai
}
write.csv(data1,file = ".../data/dnapro_f/phage_dnapro_sim_K.csv",row.names = F)

data = read.csv(".../data/species/host_gipsim.csv",header = F)
data1 = as.data.frame(data)
data1 = apply(data1, 1, as.numeric)
K = 6  ##the number of adaptive neighbors
for (i in 1:nrow(data1)){
  datai = as.numeric(data1[i,])
  datai[which(datai < sort(unique(datai),decreasing = T)[K])]=0
  data1[i,] = datai
}
write.csv(data1,file = ".../data/species/host_gipsim_K.csv",row.names = F)





################# ################# ################# ################# 
#################    ANOVA   ##################

N100 = read.csv("...data/result/N-test/N_AUC.csv",,header = T)
# 创建示例数据
group <- rep(c("A", "B", "C","D","E","F","G","H","I","J"), each = 10)
value <- c(N100$N10, N100$N20, N100$N30,N100$N40,N100$N50,N100$N60,N100$N70,N100$N80,N100$N90,N100$N100)
data <- data.frame(group, value)
result <- pairwise.t.test(data$value,data$group,p.just.method = "fdr")  #bonferroni


N310 = read.csv("...data/result/N-test/N_AUC310.csv",header = T)
# 创建示例数据
group <- rep(c("A", "B", "C","D","E","F","G","H"), each = 10)
value <- c(N310$N3, N310$N4, N310$N5,N310$N6,N310$N7,N310$N8,N310$N9,N310$N10)
data <- data.frame(group, value)
result <- pairwise.t.test(data$value,data$group,p.just.method = "fdr")  #bonferroni


anova_result <- aov(value ~ group, data = data)

print(anova_result)
summary(anova_result)