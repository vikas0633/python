### script to make hist plot

mydata = read.table("temp", header=FALSE, colClasses=c("numeric"))  # read text file

### make hist plot
hist(mydata)

