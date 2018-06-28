library(ggplot2)
library(ggsci)

data = read.table('contaminant_ratio.txt')
colnames(data) = c('id', 'len', 'con_len', 'rate')

> head(data)
           id    len con_len        rate
1 tig00024875  41173     193 0.004687538
2 tig00013505  30183       0 0.000000000
3 tig00073569  20123       0 0.000000000
4 tig00024873  57735       0 0.000000000
5 tig00000835  62928     287 0.004560768
6 tig00076896 311034    1970 0.006333713

data = within(data, {
	category <- NA
	category[rate > 0.1] <- "very_high"
	category[rate < 0.01] <- "low"
	category[rate >= 0.01 & rate <= 0.1] <- "moderate"})

> table(data$category)

      low  moderate very_high 
    18869      2075        30

ggplot(data, aes(x=len/1000, y=con_len/1000, colour=category)) + geom_point() + labs(x='Contig len(x1000)', y='Contaminated len (x1000)', colour='Contaminated Rate') + scale_color_jco(labels=c('[0, 0.01)', '[0.01, 0.1]', '[0.1, )')) + theme_classic()

save.image('contaminant_ratio.RData')
