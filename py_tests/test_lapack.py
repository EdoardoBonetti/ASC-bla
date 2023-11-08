
from ASCsoft.bla import Vector
import ASCsoft.bla as bla

# print the bla possible imports
print (dir(bla))

n = 1000000
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



