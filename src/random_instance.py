import random
import sys

n = int(sys.argv[1])
seed = int(sys.argv[2])
random.seed(seed + n * 100)
print(n)
points = []
for i in range(n):
  while True:
    x = random.randint(-100, 100)
    y = random.randint(-100, 100)
    if (x, y) not in points:
      break
  
  points.append((x, y))
  print(x,y)
