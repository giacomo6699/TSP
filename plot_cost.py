import matplotlib.pyplot as plt
import sys

y = []
file_in = open(sys.argv[1], "r")
#file_in = open("./plot/costs.txt", "r")
for x in file_in.readlines():
    if x[0:2] == "$$":
        y.append(float(x[2:-1].replace(" ", "")))
file_in.close()

plt.plot(y, marker="o")
plt.ylabel("Cost")
plt.xlabel("Iteration")
plt.title("Cost Plot")
plt.grid(visible=True)
plt.show()