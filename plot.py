import csv
from matplotlib import pyplot as plt
import sys
filename = sys.argv[1]

with open(filename) as f:
    reader = csv.reader(f)
    header = next(reader)
    num, tensor, ref = [], [], []
    for row in reader:
        num.append(row[0])
        tensor.append(row[1])
        ref.append(row[2])
plt.xlabel("test id")
plt.ylabel("Time(us)")
# plt.plot(num, tensor, label="Tensor")
plt.plot(num, ref, label="time")
plt.legend()
ax = plt.gca()
# ax.xaxis.set_major_locator(plt.MultipleLocator(1024))
ax.yaxis.set_major_locator(plt.MultipleLocator(20))
plt.savefig(filename.split(".")[0])
