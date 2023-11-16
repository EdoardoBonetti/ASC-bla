
from TomBino.bla import Vector
from TomBino.bla import Matrix
from TomBino.bla import LapackLU

# print the bla possible imports


n = 1000000
x = Vector(n)
y = Vector(n)

for i in range(n):
    x[i] = i
    y[i] = i

#print ("x =", x[0:15], "...")
#print ("y =", y[0:15], "...")

A = Matrix(3,3)
for i in [0,1,2]:
    A[i,i] = -2
    A[i,i-1] = 1
    A[i,i+1] = 1
        

z = x * y
print ("z =", z)

print( "A =" , A )

sum = 0
for i in range(n):
    sum += i*i
print ("sum =", sum)

lu = LapackLU(A)

print(lu)






