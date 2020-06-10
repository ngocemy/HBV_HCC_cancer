library(dplyr)
library(stringr)
library(ggplot2)

dir="~/Desktop/results"
files=list.files(dir)
for (f in file.path(dir,files)){
  data = read_tsv(f)
  p = ggplot(data, aes(x = "", y = observed, fill = annotation)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0)+
    theme_void()
  ggsave(p, file=paste(str_extract(f,"[:alnum:]+_1*"),".pie.png",sep = ""))
}

for (f in file.path(dir,files)){
  data = read_tsv(f)
  p_bar = ggplot(data, aes(x = annotation, y = l2fold, fill=annotation)) +
    geom_bar(stat="identity") + coord_flip() + 
    ggtitle(str_extract(f,"[:alnum:]+_1*"))
  ggsave(p_bar, file=paste(str_extract(f,"[:alnum:]+_1*"),".bar.png",sep = ""))
}

