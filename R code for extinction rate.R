metazoadata = read.csv("pbdb-dm.csv",header = TRUE)
bindata = read.table("C:\\Users\\DAI Xu\\Desktop\\Phanerozoic  metazoa diversity\\bindata2.txt",header = TRUE)
o1 = nrow(metazoadata) ##847323 occ

####bin number######
metazoadata$binno[which(metazoadata$stage == "Fortunian")] = 1
metazoadata$binno[which(metazoadata$stage == "Stage 2")] = 2
metazoadata$binno[which(metazoadata$stage == "Stage 3")] = 3
metazoadata$binno[c(which(metazoadata$stage == "Stage 4"),
                    which(metazoadata$stage == "Stage 4"))] = 4
metazoadata$binno[c(which(metazoadata$stage == "Drumian"),
                    which(metazoadata$stage == "Guzhangian"),
                    which(metazoadata$stage == "Paibian"))] = 5
metazoadata$binno[c(which(metazoadata$stage == "Jiangshanian"),
                    which(metazoadata$stage == "Stage 10"))] = 6
metazoadata$binno[which(metazoadata$stage == "Tremadocian")] = 7
metazoadata$binno[c(which(metazoadata$stage == "Floian"),
                    which(metazoadata$stage == "Dapingian"))] = 8
metazoadata$binno[c(which(metazoadata$stage == "Darriwilian"),
                    which(metazoadata$stage == "Sandbian"))] = 9
metazoadata$binno[c(which(metazoadata$stage == "Katian"),
                    which(metazoadata$stage == "Hirnantian"))] = 10
metazoadata$binno[c(which(metazoadata$stage == "Rhuddanian"),
                    which(metazoadata$stage == "Aeronian"),
                    which(metazoadata$stage == "Telychian"))] = 11
metazoadata$binno[c(which(metazoadata$stage == "Sheinwoodian"),
                    which(metazoadata$stage == "Homerian"))] = 12
metazoadata$binno[c(which(metazoadata$stage == "Pridoli"),
                  which(metazoadata$stage == "Gorstian"),
                  which(metazoadata$stage == "Ludfordian"))] = 13
metazoadata$binno[c(which(metazoadata$stage == "Lochkovian"),
                    which(metazoadata$stage == "Pragian"))] = 14
metazoadata$binno[which(metazoadata$stage == "Emsian")] = 15
metazoadata$binno[c(which(metazoadata$stage == "Eifelian"),
                    which(metazoadata$stage == "Givetian"))] = 16
metazoadata$binno[which(metazoadata$stage == "Frasnian")] = 17
metazoadata$binno[which(metazoadata$stage == "Famennian")] = 18
metazoadata$binno[which(metazoadata$stage == "Tournaisian")] = 19
metazoadata$binno[(which(metazoadata$stage == "Visean") && (metazoadata$earlyage+metazoadata$earlyage)/2 > 338.8)] = 20
metazoadata$binno[(which(metazoadata$stage == "Visean") && (metazoadata$earlyage+metazoadata$earlyage)/2 < 338.8)] = 21
metazoadata$binno[which(metazoadata$stage == "Serpukhovian")] = 21
metazoadata$binno[which(metazoadata$stage == "Bashkirian")] = 22
metazoadata$binno[c(which(metazoadata$stage == "Moscovian"),
                    which(metazoadata$stage == "Kasimovian"))] = 23
metazoadata$binno[which(metazoadata$stage == "Gzhelian")] = 24
metazoadata$binno[c(which(metazoadata$stage == "Asselian"),
                    which(metazoadata$stage == "Sakmarian"))] = 25
metazoadata$binno[which(metazoadata$stage == "Artinskian")] = 26
metazoadata$binno[c(which(metazoadata$stage == "Kungurian"),
                    which(metazoadata$stage == "Roadian"))] = 27
metazoadata$binno[c(which(metazoadata$stage == "Wordian"),
                    which(metazoadata$stage == "Capitanian"))] = 28
metazoadata$binno[c(which(metazoadata$stage == "Wuchiapingian"),
                    which(metazoadata$stage == "Changhsingian"))] = 29
metazoadata$binno[c(which(metazoadata$stage == "Induan"),
                    which(metazoadata$stage == "Olenekian"))] = 30
metazoadata$binno[c(which(metazoadata$stage == "Anisian"),
                    which(metazoadata$stage == "Ladinian"))] = 31
metazoadata$binno[which(metazoadata$stage == "Carnian")] = 32
metazoadata$binno[which(metazoadata$stage == "Norian")] = 33
metazoadata$binno[which(metazoadata$stage == "Rhaetian")] = 34
metazoadata$binno[c(which(metazoadata$stage == "Hettangian"),
                  which(metazoadata$stage == "Sinemurian"))] = 35
metazoadata$binno[which(metazoadata$stage == "Pliensbachian")] = 36
metazoadata$binno[c(which(metazoadata$stage == "Toarcian"),
                    which(metazoadata$stage == "Aalenian"))] = 37
metazoadata$binno[c(which(metazoadata$stage == "Bajocian"),
                    which(metazoadata$stage == "Bathonian"),
                    which(metazoadata$stage == "Callovian"))] = 38
metazoadata$binno[which(metazoadata$stage == "Oxfordian")] = 39
metazoadata$binno[c(which(metazoadata$stage == "Kimmeridgian"),
                    which(metazoadata$stage == "Tithonian"))] = 40
metazoadata$binno[c(which(metazoadata$stage == "Berriasian"),
                    which(metazoadata$stage == "Valanginian"))] = 41
metazoadata$binno[c(which(metazoadata$stage == "Hauterivian"),
                    which(metazoadata$stage == "Barremian"))] = 42
metazoadata$binno[which(metazoadata$stage == "Aptian")] = 43
metazoadata$binno[which(metazoadata$stage == "Albian")] = 44
metazoadata$binno[which(metazoadata$stage == "Cenomanian")] = 45
metazoadata$binno[c(which(metazoadata$stage == "Turonian"),
                    which(metazoadata$stage == "Coniacian"),
                    which(metazoadata$stage == "Santonian"))] = 46
metazoadata$binno[which(metazoadata$stage == "Campanian")] = 47
metazoadata$binno[which(metazoadata$stage == "Maastrichtian")] = 48
metazoadata$binno[c(which(metazoadata$stage == "Danian"),
                    which(metazoadata$stage == "Selandian"),
                    which(metazoadata$stage == "Thanetian"))] = 49
metazoadata$binno[which(metazoadata$stage == "Ypresian")] = 50
metazoadata$binno[c(which(metazoadata$stage == "Lutetian"))] = 51
metazoadata$binno[c(which(metazoadata$stage == "Bartonian"),
                    which(metazoadata$stage == "Priabonian"))] = 52
metazoadata$binno[c(which(metazoadata$stage == "Rupelian"),
                    which(metazoadata$stage == "Chattian"))] = 53
metazoadata$binno[c(which(metazoadata$stage == "Aquitanian"),
                    which(metazoadata$stage == "Burdigalian"))] = 54
metazoadata$binno[c(which(metazoadata$stage == "Langhian"),
                    which(metazoadata$stage == "Serravallian"),
                    which(metazoadata$stage == "Tortonian"),
                    which(metazoadata$stage == "Messinian"))] = 55
metazoadata$binno[c(which(metazoadata$stage == "Zanclean"),
                    which(metazoadata$stage == "Piacenzian"))] = 56

metazoadata1 = metazoadata[-c(which(metazoadata$class == "Ostracoda"),
                             which(metazoadata$class == "Arachnida"),
                             which(metazoadata$class == "Insecta"),
                             which(metazoadata$class == "Reptilia"),
                             which(metazoadata$genus == "")),]

library(divDyn)
result = divDyn(metazoadata1,tax = "genus", bin = "binno")
write.csv(result, file = "extinction-result.csv")


###########################################
###### GF extinction rate estimation ######
###########################################

div_data = read.csv("extinction-result.csv.csv",header = TRUE) #load data
est_gf = function(k){
  a = (div_data$t2d[k]+div_data$tPart[k])
  b = (div_data$t3[k]+div_data$tPart[k]+div_data$tGFu[k])
  c = a+b
  erate = rep(NA, 10000)
  for (i in 1:10000){
    sam_data = sample(c, size = c, replace = TRUE)
    erate[i] = log(length(which(sam_data<=a))/length(which(sam_data > a)))
  }
  return(c(quantile(erate,0.025),quantile(erate,0.05),
           quantile(erate,0.15),mean(erate),
           quantile(erate,0.85),quantile(erate,0.95),
           quantile(erate,0.975),sd(erate)))
}
gf_erate = matrix(NA, 56, 8)
for (i in 3:55){
  gf_erate[i,] = est_gf(i)
  print(i)
}
write.csv(gf_erate, file = "gf_erate.csv")

