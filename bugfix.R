bla = sort(as.character(toptrans$transition_group_id))
bla[1:10]
write.table(bla,file="xx.txt")
bla = read.table("xx.txt",stringsAsFactors=F)
bla = bla[,1]
bla[1:10]
sort(bla)[1:10]
