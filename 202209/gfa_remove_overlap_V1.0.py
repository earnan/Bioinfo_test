#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   gfa_remove_overlap
#         Author:   yujie
#    Description:   gfa_remove_overlap
#        Version:   1.0
#           Time:   2022/08/13 18:58:15
#  Last Modified:   2022/08/13 18:58:15
#        Contact:   hi@arcsona.cn
#        License:   Copyright (C) 2022
#
##########################################################
from Bio import SeqIO
from Bio.Seq import Seq
from icecream import ic
import argparse
import linecache
import os
import pretty_errors
import re
import sys
import time
import copy

parser = argparse.ArgumentParser(
    add_help=False, usage='\
\npython3   gfa_remove_overlap_V1.0.py\n\
功能：\n\
\n\
1.常规使用\n\
1.1 -i [ ] -o [ ] \n\
\n\
2.其他使用\n\
2.1 -i [ ] -o [ ] \n\
\n\
Path: E:\OneDrive\jshy信息部\Script\Bioinfo_test\gfa_remove_overlap_V1.0.py\n\
Path: /share/nas1/yuj/script/Bioinfo_test/gfa_remove_overlap_V1.0.py\n\
Version: 1.0'
)
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--infile', metavar='[infile]', help='infile', type=str, default='E:\\OneDrive\\jshy信息部\\Script\\Bioinfo_test\\001_best_spades_graph.gfa', required=False)
optional.add_argument(
    '-o', '--outfile', metavar='[outfile]', help='outfile', type=str, default='F:/', required=False)
optional.add_argument('-c1', '--flag1', help='run step 1?默认是,不运行则-c1',
                      action='store_false', required=False)
optional.add_argument('-c2', '--flag2', help='run step 2?默认否,运行则-c2 ',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()


# ####################################################################
# 读取文件
ingfa_path = args.infile
with open(ingfa_path, 'r') as ingfa_handle:
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
                # 39 + 7:+,8:+    对应+链末尾127bp
                l_line_dict[index1][0].append(index2+':'+index2_strand)
            elif index1_strand == '-':
                l_line_dict[index1][1].append(
                    index2+':'+index2_strand)  # 39 - 4:+  对应+链开头127bp
l_line_dict = dict(
    sorted(l_line_dict.items(), key=lambda x: x[0], reverse=False))  # 排序

l_line_dict_1 = copy.deepcopy(l_line_dict)
ic('未去重', l_line_dict_1)

# ###########################################################
# 移除重叠碱基 子函数


def remove_seq_overlap(index_strand_re_pos, seq, ovl):  # 移除重叠碱基
    if index_strand_re_pos == 'end':  # 移除末尾
        new_seq = seq[:-ovl]
    elif index_strand_re_pos == 'start':
        new_seq = seq[ovl:]  # 移除开头
    return new_seq

# ############################################################################
# 获取关联的overlap关系 子函数


def append_new_ovl_list(i, k, long_seq_id, index1_strand_re_pos, new_ovl_list):
    index2 = int(re.findall(r'\d+', i)[0])
    index2_strand = i.split(':')[-1]
    if index2_strand == '+':
        index2_strand_re_pos = 'start'
    else:
        index2_strand_re_pos = 'end'
    if index2 == long_seq_id:
        index1 = k
        new_ovl_list.append([
            str(index1)+':'+index1_strand_re_pos, str(index2)+':'+index2_strand_re_pos, i])
    return 0


def get_long_seq_new_ovl_list(long_seq_id, l_line_dict):
    index1_strand_re_pos = ''
    new_ovl_list = []
    for k, v in l_line_dict.items():
        if len(v[0]) > 0:
            index1_strand_re_pos = 'end'
            for i in v[0]:
                append_new_ovl_list(i, k, long_seq_id,
                                    index1_strand_re_pos, new_ovl_list)
        if len(v[1]) > 0:
            index1_strand_re_pos = 'start'
            for i in v[1]:
                append_new_ovl_list(i, k, long_seq_id,
                                    index1_strand_re_pos, new_ovl_list)
    return new_ovl_list  # 数值,字符,数值,字符开头末尾

# ###################################################################
# 移除关联的overlap


def remove_new_ovl(new_ovl_list, long_seq_id_re_pos, s_line_dict, l_line_dict):
    if len(new_ovl_list) > 0:
        for tmp_i in new_ovl_list:
            new_index1 = int(re.findall(r'\d+', tmp_i[0])[0])
            new_index1_strand_re_pos = tmp_i[0].split(':')[-1]
            new_index2 = int(re.findall(r'\d+', tmp_i[1])[0])
            new_index2_strand_re_pos = tmp_i[1].split(':')[-1]
            i_content = tmp_i[2]
            if new_index2_strand_re_pos != long_seq_id_re_pos:
                if len(s_line_dict[new_index1]) > len(s_line_dict[new_index2]):
                    long_seq_id = new_index1
                    long_seq_id_re_pos = new_index1_strand_re_pos
                else:
                    long_seq_id = new_index2
                    long_seq_id_re_pos = new_index2_strand_re_pos

                s_line_dict[long_seq_id] = remove_seq_overlap(
                    long_seq_id_re_pos, s_line_dict[long_seq_id], ovl)  # 移除序列里的重叠碱基
                if new_index1_strand_re_pos == 'end':
                    l_line_dict[new_index1][0].remove(
                        i_content)  # 移除字典里的重叠关系
                elif new_index1_strand_re_pos == 'start':
                    l_line_dict[new_index1][1].remove(
                        i_content)  # 移除字典里的重叠关系
    return 0

# ##################################################################
# 按顺序,移除出现的overlap


def remove_edit(tmp_list, s_line_dict, index1, index1_strand_re_pos, l_line_dict):
    for i in tmp_list:
        index2 = int(i.split(':')[0])  # 数值
        index2_strand = i.split(':')[1]  # 方向+-
        if len(s_line_dict[index1]) > len(s_line_dict[index2]):
            long_seq_id = index1
            long_seq_id_re_pos = index1_strand_re_pos
        else:
            long_seq_id = index2
            if index2_strand == '+':
                long_seq_id_re_pos = 'start'
            else:
                long_seq_id_re_pos = 'end'
        s_line_dict[long_seq_id] = remove_seq_overlap(
            long_seq_id_re_pos, s_line_dict[long_seq_id], ovl)  # 移除序列里的重叠碱基
        tmp_list.remove(i)  # 移除字典里的重叠关系
        '''继续删除相关序列'''
        new_ovl_list = get_long_seq_new_ovl_list(
            long_seq_id, l_line_dict)  # 根据上一个长id寻找其与其他序列的重叠关系
        remove_new_ovl(new_ovl_list, long_seq_id_re_pos,
                       s_line_dict, l_line_dict)
    return 0


# ##############################################################################################################
# 思路2  直接删除ovl    留下交汇结点  其他的删掉

for index1, v in l_line_dict.items():
    '''list1相关 index1 +    index2 +/- '''
    if len(v[0]) > 0:  # index1 +    index2 +/-
        tmp_list = v[0]
        index1_strand_re_pos = 'end'
        remove_edit(tmp_list, s_line_dict, index1,
                    index1_strand_re_pos, l_line_dict)
    '''list2相关 index1 -    index2 +/- '''
    if len(v[1]) > 0:
        tmp_list = v[1]
        index1_strand_re_pos = 'start'
        remove_edit(tmp_list, s_line_dict, index1,
                    index1_strand_re_pos, l_line_dict)

ic('去重', l_line_dict)

# ##############################################################################
# 输出新的gfa
outgfa_path = "E:\\OneDrive\\jshy信息部\\Script\\Bioinfo_test\\edit-001_best_spades_graph.gfa"
with open(outgfa_path, 'w') as outgfa:
    for k, v in s_line_dict.items():
        outgfa.write('S\t{}\t{}\n'.format(k, v))
    print(l_line_dict_1)
    for k, v in l_line_dict_1.items():
        if len(v[0]) > 0:
            for i in v[0]:
                outgfa.write(
                    'L\t{}\t+\t{}\t{}\t0M\n'.format(k, i.split(':')[0], i.split(':')[-1]))
        if len(v[1]) > 0:
            for i in v[1]:
                outgfa.write(
                    'L\t{}\t-\t{}\t{}\t0M\n'.format(k, i.split(':')[0], i.split(':')[-1]))
