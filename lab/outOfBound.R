
rm (list = ls())
myArray = 1:2
myArray [0]
myArray [-1]
myIndex = -1
myArray [myIndex]
myArray [3]

myDoubleArray = array (data = 1:3, dim = c(2, 4))
myDoubleArray
myDoubleArray [20]
myDoubleArray [0]
myDoubleArray [c(0)]
myDoubleArray [c(1)]
myDoubleArray [c(1, 2)]
class (NULL)
class (NA)
# only here do we have error
myDoubleArray [20, 20]

