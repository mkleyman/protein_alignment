library(ggplot2)
library(dplyr)
library(ggrepel)
colnames(search_table) = c("curve","score","threshold", "mouse", "human")
search_table_fixed = mutate(search_table, human_years = human/365.0)
rate_thresh = filter(search_table_fixed, threshold == 0.8)
table(search_table$curve)

ggplot(rate_thresh, aes(x=mouse,y=human_years, color=curve))+
  geom_line(size=2)+xlab("Mouse(days) since conception")+ylab("Human(years) since birth")+
  ylim(0,10)+ggtitle("Pearson Correlation Threshold 0.8")+theme(text = element_text(size=26))


colnames(go_table_uniq) = c("go_category","significant","total")
table(go_table_uniq$significant)
go_table_significant = filter(go_table_uniq, significant>5)
ggplot(go_table_significant, aes(x=total,y=significant, label=go_category))+geom_point(size=3)+
  geom_text_repel(size=6)+ggtitle("Significant GO categories for Quadratic Function at .7 Threshold")+
  theme(text = element_text(size=26))

go_table_exp_uniq <- read.csv("~/protein_alignment/processed/go_table_exp_uniq.csv", header=FALSE, stringsAsFactors=FALSE)

colnames(go_table_exp_uniq) = c("go_category","significant","total")
go_table_exp_significant = filter(go_table_exp_uniq, significant>5)
ggplot(go_table_exp_significant, aes(x=total,y=significant, label=go_category))+geom_point(size=3)+
  geom_text_repel(size=6)+ggtitle("Significant GO categories for Exponential Function at .7 Threshold")+
  theme(text = element_text(size=26))

colnames(significance_results) = c("data_source","found_homologs")
ggplot(significance_results, aes(x=found_homologs,fill=data_source))+geom_histogram(size=2)+
  theme(text = element_text(size=26))+ggtitle("Number of Homologs found in Random Selection of Proteins")


chosen_proteins <- read.csv("~/protein_alignment/processed/chosen_proteins.txt", header=FALSE, stringsAsFactors=FALSE)
colnames(chosen_proteins) = c("mouse","human","type","time","expression")
table(chosen_proteins$mouse)

homolog = filter(chosen_proteins, mouse=="Q9Z1N5")
homolog[1,]
mouse_table = filter(homolog, type=="mouse")
ggplot(mouse_table, aes(x=time,y=expression))+geom_point(size=4)+geom_smooth()+
  theme(text = element_text(size=26))+ggtitle("Q9Z1N5 (Mouse)")

human_table = filter(homolog, type=="human")
ggplot(human_table, aes(x=time,y=expression))+geom_point(size=4)+geom_smooth()+
  theme(text = element_text(size=26))+ggtitle("Q13838 (Human)")

align_table = filter(homolog, type!="mouse" & type!="human" )
ggplot(align_table, aes(x=time,y=expression, color=type))+geom_point(size=4)+geom_smooth(size=1)+
  theme(text = element_text(size=26))+ggtitle(" Q9Z1N5 (Mouse) and Q13838 (Human)")

colnames(exp) = c("score","thresh","co1","co2")
summarise(group_by(exp,thresh),mean_score = mean(score, na.rm = TRUE))

linear <- read.csv("~/protein_alignment/processed/FDR2/linear.csv", header=FALSE, stringsAsFactors=FALSE)
colnames(linear) = c("score","thresh","co1","co2")
summarise(group_by(linear,thresh),mean_score = mean(score, na.rm = TRUE))


sqrt <- read.csv("~/protein_alignment/processed/FDR2/sqrt.csv", header=FALSE, stringsAsFactors=FALSE)
colnames(sqrt) = c("score","thresh","co1","co2")
summarise(group_by(sqrt,thresh),mean_score = mean(score, na.rm = TRUE))

lgthm <- read.csv("~/protein_alignment/processed/FDR2/logarithm.csv", header=FALSE, stringsAsFactors=FALSE)
colnames(lgthm) = c("score","thresh","co1","co2")
summarise(group_by(lgthm,thresh),mean_score = mean(score, na.rm = TRUE))


diff_proteins <- read.csv("~/protein_alignment/processed/diff_proteins.txt", header=FALSE, stringsAsFactors=FALSE)
colnames(diff_proteins) = c("mouse","human","type","time","expression")
table(diff_proteins$mouse)
good = c("Q6PGH2", "Q91XF0", "P14824")
bad = c("P20152","Q91XF0 ")
homolog = filter(diff_proteins, mouse=="P20152")
homolog[1,]
mouse_table = filter(homolog, type=="mouse")
ggplot(mouse_table, aes(x=time,y=expression))+geom_point(size=4)+geom_smooth()+
  theme(text = element_text(size=26))+ggtitle("Q91XF0 (Mouse)")

human_table = filter(homolog, type=="human")
ggplot(human_table, aes(x=time,y=expression))+geom_point(size=4)+geom_smooth()+
  theme(text = element_text(size=26))+ggtitle("Q9H910(Human)")

align_table = filter(homolog, type!="mouse" & type!="human" )
ggplot(align_table, aes(x=time,y=expression, color=type))+geom_point(size=4)+geom_smooth(size=1)+
  theme(text = element_text(size=26))+ggtitle(" Q91XF0 (Mouse) and Q9H910 (Human)")

xvals = c(1,2,3,4,5,6,7,8,9,10)
yvals= c(3,5,6,9,10,11,13,-5,0,8)
spl = smooth.spline(xvals,yvals, all.knots = FALSE, nknots = 4)
spl$lambda