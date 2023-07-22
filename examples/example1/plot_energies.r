setwd("~/Docs/tetmag/example1/")
library(ggplot2)
logdata <- read.table("disk.log")
p <- ggplot(logdata,aes(x=V1)) + geom_line(aes(y=V2, color = "V2"))+
  geom_line(aes(y=V3, color = "V3"))+ geom_line(aes(y=V4, color = "V4"))+
  scale_color_manual(values = c("red", "darkgreen", "blue"), labels=c("total", "demag", "exchange")) +
  scale_y_log10() + xlab("time [ps]") + ylab(expression(paste ("energy density [", J/m^{3}, "]")))+ theme(legend.title = element_blank()) 
print(p)