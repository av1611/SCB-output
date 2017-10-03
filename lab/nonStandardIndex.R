

myMatrix = matrix (data = 1:9, ncol = 3, nrow = 3)
myMatrix
myMatrix$`1`
myDataframe = as.data.frame(myMatrix)
myDataframe
myDataframe$V1
# myDataframe [1, 1]
# myDataframe [1, ]
# myDataframe [, 1]
colnames (myDataframe) = c (1, 2, 3)
myDataframe$`1`
myIndex = 2
myExpr <- paste0("myDataframe$`", myIndex, "`")
myDataframe
eval(parse(text = myExpr))




