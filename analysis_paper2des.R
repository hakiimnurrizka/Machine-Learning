### Analysis of student's intention to participate in online class 
library(readxl)
library(ggplot2)
library(nortest)
library(car)
data_paper2des <- read_excel("G:/My Drive/Github/Machine-Learning/data.paper2des.xlsx")

ggplot(data = data_paper2des, aes(x = int.score, y = posting)) + 
  geom_point() + plot(data_paper2des$int.score,data_paper2des$posting) +
  xlab("Intention score") + ylab("Number of posting")

shapiro.test(data_paper2des$int.score)
shapiro.test(data_paper2des$posting)
shapiro.test(data_paper2des$grade)

cor.test(data_paper2des$int.score, data_paper2des$posting, method = "kendall")
cor.test(data_paper2des$int.score, data_paper2des$posting, method = "spearman")

ggplot(data = data_paper2des, aes(x = int.score, y = grade)) + 
  geom_point() + plot(data_paper2des$int.score,data_paper2des$grade) +
  ylab("Final grade") + xlab("Intention score")
cor.test(data_paper2des$int.score, data_paper2des$grade, method = "kendall")
cor.test(data_paper2des$int.score, data_paper2des$grade, method = "spearman")

ggplot(data = data_paper2des, aes(x = posting, y = grade)) + 
  geom_point() + plot(data_paper2des$posting,data_paper2des$grade) +
  ylab("Final grade") + xlab("Number of posting")
cor.test(data_paper2des$grade, data_paper2des$posting, method = "kendall")
cor.test(data_paper2des$grade, data_paper2des$posting, method = "spearman")


lm.2 = lm(grade ~ poly(posting,2, raw =T), data_paper2des)
summary(lm.2)
lm.3 = lm((grade) ~ -1+ (sin(log(int.score))) + ((posting)), data_paper2des[,])

summary(lm.3)
par(mfrow = c(1,1))
plot(lm.3)
vif(lm.3)
