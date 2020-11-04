import subprocess
import csv
import sys
filename = sys.argv[1]

# print("_data/"+filename+"_beg_pos.bin")  com-orkut.txt_b

beg_pos="_data/"+filename+"_beg_pos.bin"
csr="_data/"+filename+"_csr.bin"
weight="_data/"+filename+"_weight.bin"
# print(weight)

def get_time_num(num):
    # out = subprocess.check_output(
    #     ['./main.out', str(num), filename, '2']).decode('utf-8')
    out = subprocess.check_output(
        ['./graph_partition.out', beg_pos, csr, weight]).decode('utf-8')
    temp = out.split("\n")
    print(temp)

    time = float(temp[-1])
    edge_after = int(temp[-2])
    edge_before = int(temp[-3])
    
    return time, edge_after, edge_before

# get_time_num(1)

with open(f'{filename}_num_1_3_1.csv', 'w+') as f:
    headers = ['time', 'edge_after', 'edge_before']
    f_csv = csv.writer(f)
    f_csv.writerow(headers)
    for i in range(1, 10, 1):
        # if():
        #     print(i)
        f_csv.writerow(get_time_num(i))
