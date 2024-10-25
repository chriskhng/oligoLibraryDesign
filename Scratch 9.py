import math

pattern = 's..s..s'
inputPlusMinus = 15
spanPlusMinus = inputPlusMinus - math.ceil(len(pattern)/2)+1
print(spanPlusMinus)


