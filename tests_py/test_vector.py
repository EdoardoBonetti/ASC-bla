import TomBino as tb
from TomBino.bla import Vector, VectorInt , VectorComplex

print(dir(tb.bla))



from TomBino.bla import *

n = 100
x = Vector(n)
y = Vector(n)

for i in range(n):
    x[i] = i
    y[i] = i

#print ("x =", x[0:15], "...")
#print ("y =", y[0:15], "...")


z = x * y
print ("z =", z)

sum = 0
for i in range(n):
    sum += i*i
print ("sum =", sum)

## do the same for int and complex
x_Int = VectorInt(n)
y_Int = VectorInt(n)

x_Complex = VectorComplex(n)
y_Complex = VectorComplex(n)

for i in range(n):
    x_Int[i] = i
    y_Int[i] = i

    x_Complex[i] = i
    y_Complex[i] = i
#print ("x =", x[0:15], "...")
#print ("y =", y[0:15], "...")
#print ("x_Int =", x_Int[0:15], "...")
#print ("y_Int =", y_Int[0:15], "...")

z_Int = x_Int * y_Int
print ("z_Int =", z_Int)

z_Complex = x_Complex * y_Complex
print ("z_Complex =", z_Complex)






