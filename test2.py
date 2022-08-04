from Bio import SeqIO
from Bio.Seq import Seq
from icecream import ic
import argparse
import linecache
import os
import re
import time
import sys
parser = argparse.ArgumentParser(
    add_help=False, usage='\
\npython3   test.py\n\
step1\n\
step2\n\
V1.0')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--infile', metavar='[infile]', help='infile', type=str, default='E:/', required=False)
optional.add_argument(
    '-o', '--outfile', metavar='[outfile]', help='outfile', type=str, default='F:/', required=False)
optional.add_argument('-c1', '--flag1', help='run step 1?默认是,不运行则-c1',
                      action='store_false', required=False)
optional.add_argument('-c2', '--flag2', help='run step 2?默认否,运行则-c2 ',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()

# ###########################################################
# 子函数


def remove_overlap(index1_strand, seq, ovl):  # 移除重叠部分
    if index1_strand == '+':  # 移除末尾
        new_seq = seq[:-ovl]
    elif index1_strand == '-':
        new_seq = seq[ovl:]  # 移除开头
    return new_seq


# ####################################################################
# 读取文件
ingfa_path = "E:\\OneDrive\\jshy信息部\\Bioinfo_test\\001_best_spades_graph.gfa"
with open(ingfa_path, 'r') as ingfa_handle:
    '''定义变量'''
    s_line_dict = {}
    l_line_dict = {}
    num_list = []

    for line in ingfa_handle:
        '''id和对应序列都存入字典s'''
        if line.startswith('S'):
            s_content = line.strip().split()
            s_line_dict[int(s_content[1])] = s_content[2]

            '''overlap关系存入字典l'''
        elif line.startswith('L'):
            content = line.strip().split()
            index1 = int(content[1])  # 1链id
            index1_strand = content[2]  # 1链方向
            index2 = content[3]  # 2链id
            index2_strand = content[4]  # 2链方向
            ovl = int(re.findall(r'\d+', content[5])[0])  # 127
            if index1 not in num_list:
                l_line_dict[index1] = [[], []]  # + -
                num_list.append(index1)

            '''存储1链    与  对应的2链及2链方向'''
            if index1_strand == '+':
                l_line_dict[index1][0].append(index2+':'+index2_strand)
                # list1 = []  # 39 + 7:+,8:+    对应+链末尾127bp
            elif index1_strand == '-':
                l_line_dict[index1][1].append(index2+':'+index2_strand)
                # list2 = []  # 39 - 4:+  对应+链开头127bp
ic(l_line_dict)
# #######################################################################

# s_line_dict  序列
# l_line_dict   重叠关系

# 初始化
double_ovl_index1 = []  # 双端ovl
single_ovl_index1 = []  # 单端ovl
l_line_dict = dict(
    sorted(l_line_dict.items(), key=lambda x: x[0], reverse=False))  # 排序
for k, v in l_line_dict.items():
    if len(v[0]) > 0 and len(v[1]) > 0:
        double_ovl_index1.append(k)
    else:
        single_ovl_index1.append(k)

# 去重
for i in single_ovl_index1:
    for m in l_line_dict[i][0]:  # k = 39:+
        if int(re.findall(r'\d+', m)[0]) in double_ovl_index1:
            l_line_dict[i][0].remove(m)
    for n in l_line_dict[i][1]:  # v = 3:+
        if int(re.findall(r'\d+', n)[0]) in double_ovl_index1:
            l_line_dict[i][1].remove(n)

# ic(l_line_dict)
# ######################################################################
# 修改序列
# s_line_dict  旧的序列
# l_line_dict   新的重叠关系

# for k, v in s_line_dict.items():
#print(k, len(v))

double_new_ovl_index1 = []  # 双端ovl
single_new_ovl_index1 = []  # 单端ovl
for k, v in l_line_dict.items():
    if len(v[0]) > 0:
        s_line_dict[k] = remove_overlap('+', s_line_dict[k], ovl)
    if len(v[1]) > 0:
        s_line_dict[k] = remove_overlap('-', s_line_dict[k], ovl)

print('-----------------------------------------')
# for k, v in s_line_dict.items():
#print(k, len(v))
