import curves

print "Quadratic:\n\n"
#quadratic best
coefs = [-1390.99,213,-2.42]

for i in xrange(16,44):
    print curves.quadratic(i,coefs)/365.0

print "Cubic:\n\n"

#cubic best
coefs =[-4.45953766e+03,   3.72329374e+02,  -4.40245508e+00,   3.76569485e+03]
for i in xrange(16,44):
    print curves.quadratic(i,coefs)/365.0