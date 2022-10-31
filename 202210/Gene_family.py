#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   Gene_family.py
#         Author:   yujie
#    Description:   Gene_family.py
#        Version:   1.0
#           Time:   2022/10/29 15:36:11
#  Last Modified:   2022/10/31 15:36:11
#        Contact:   hi@arcsona.cn
#        License:   GNU General Public License v3.0
#
##########################################################
#from Bio import SeqIO
#from Bio.Seq import Seq
# from humre import *  # 正则
# from icecream import ic  # 打印
import argparse  # 命令行
import linecache  # 大文件行读取
import os  # 目录路径
# import pretty_errors  # 错误提示
import re  # 正则
import sys
import time
import copy  # 深度拷贝
#import pandas as pd
#import numpy as np
#import matplotlib.pyplot as plt
parser = argparse.ArgumentParser(
    add_help=False, usage='\n\
\n\
##########################################################\n\
#\n\
#       Filename:   Gene_family.py\n\
#         Author:   yujie\n\
#    Description:   Gene_family.py\n\
#        Version:   1.0\n\
#           Time:   2022/10/11 10:49:39\n\
#  Last Modified:   2022/10/31 15:40:04\n\
#        Contact:   hi@arcsona.cn\n\
#        License:   GNU General Public License v3.0\n\
#\n\
##########################################################\n\
\n\
\npython3   Gene_family.py\n\
Function:\n\
1.常规使用\n\
1.1 -i [resultfile.cafe] -o [stat.xls] \n\
\n\
##########################################################\n\
Path: E:\OneDrive\jshy信息部\Script\Bioinfo_test\202210\Gene_family.py\n\
Path: /share/nas1/yuj/script/Bioinfo_test/Gene_family.py\n\
Version: 1.0\n\
##########################################################\n\
'
)
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--innwk', metavar='[nwk]', help='输入文件，默认为该程序路径下的文件', type=str, default="E:\\OneDrive\\jshy信息部\\Script\\Bioinfo_test\\202210\\resultfile.cafe", required=False)
optional.add_argument(
    '-o', '--outtab', metavar='[table]', help='输出文件，默认为该程序路径下的文件', type=str, default="E:\\OneDrive\\jshy信息部\\Script\\Bioinfo_test\\202210\\statistics.xls", required=False)
#optional.add_argument('-c1', '--flag1', help='run step 1?默认是,不运行则-c1',action='store_false', required=False)
#optional.add_argument('-c2', '--flag2', help='run step 2?默认否,运行则-c2 ',action='store_true', required=False)
optional.add_argument('-info', help='更新日志,使用时-info',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help', help='[help_information]')
args = parser.parse_args()

if args.info:
    print('\n更新日志:')
    print('\t2022/10/11  首次创建文件')
    print('\t2022/10/29  修改文件写入的内容')
    print('\t2022/10/31  优化一些可能出bug的代码，增加完整的信息提示')
    sys.exit(0)

# 读取输入文件，获取最终分支id  调用一次即可


def read_nwk(innwk_path):
    with open(innwk_path, 'r') as innwk_handle:
        for line in innwk_handle:
            if line.startswith('# IDs of nodes:'):
                ori_list_nodes_ids = re.findall(
                    r'[A-Za-z0-9<>]+', line.split(':')[1])  # 取出 字母数字<>
                filtered_list_nodes_ids = []
                for i in ori_list_nodes_ids:
                    if not i.startswith('<'):
                        filtered_list_nodes_ids.append(i)
                simplified_list_nodes_ids = re.findall(
                    r'[A-Za-z]+', line.split(':')[1])  # 取出字母
                break
    return simplified_list_nodes_ids, filtered_list_nodes_ids  # 一个为li  一个为li<0>

# 传进去一行  需要每次调用


def get_list_node(ori_line):
    # \w 等价于[A-Za-z0-9_] 字母数字下划线
    # li_2', '32', '3768', 'taozi_0', '7', '686', 'Pmu_5', '7', '686', '_2'
    ori_list_node = re.findall(r'[\w]+', ori_line)
    n = 0
    list_node = []
    for i in ori_list_node:
        if i.find('_') >= 0:  # 找出带有数量 实际上有意义的字段（节点）
            n += 1
            list_node.append(str(n)+','+i)  # 会有重复字段，加上序号就能分开了
    return list_node  # 1,li_36


if __name__ == '__main__':
    #################################################################
    # 格式化成2016-03-20 11: 45: 39形式
    start_time = time.time()
    print('Start Time : {}'.format(time.strftime(
        '%Y-%m-%d %H:%M:%S', time.localtime())))
    #################################################################

    print('\n---------------------------Step 00 initialization------------------------------------')
    # 文件初始化
    innwk_path = args.innwk
    outtab_path = args.outtab
    # 动态创建变量   函数初始化
    createVar = locals()
    print('Done')

    print('\n---------------------------Step 01 get_the_final_node_id_and_dynamically_generate_the_required_dictionary------------------------------------')
    simplified_list_nodes_ids, filtered_list_nodes_ids = read_nwk(innwk_path)
    print(filtered_list_nodes_ids)
    print(simplified_list_nodes_ids)
    for simplified_node_id in simplified_list_nodes_ids:
        createVar['dict_'+simplified_node_id+'_allinfo'] = {}  # 该物种字典,用于写入
        createVar['dict_'+simplified_node_id+'_allinfo']['shrunk'] = []
        createVar['dict_'+simplified_node_id+'_allinfo']['expanded'] = []
        createVar['dict_'+simplified_node_id+'_allinfo']['unchanged'] = []
    print('Done')

    print('\n---------------------------Step 02 read_the_file_line_by_line_for_extraction------------------------------------')
    innwk_contents = linecache.getlines(innwk_path)
    skip_row_count = 0
    for line in innwk_contents:
        skip_row_count += 1
        if line.startswith("'ID'"):
            break
    print("skip the first {} lines".format(skip_row_count))
    innwk_contents = innwk_contents[skip_row_count:]  # 此例为跳过前十行   从有进化树的一行开始
    valid_row_count = 0
    trigger_flag = False
    for line in innwk_contents:
        if line.strip() == '':
            print('{}_rows_in_total_are_valid'.format(
                valid_row_count))  # 此if条件在使用linecache时不会触发  linecache只会读取非空行  open()读取方式会触发
            trigger_flag = True
        else:
            valid_row_count += 1
            family_name = line.split('\t')[0]
            ori_line = line.split('\t')[1]
            list_node = get_list_node(ori_line)
            # 获取node在该行的index   需要每次调用
            ori_filtered_node_index = 0  # node在该行的绝对索引
            ori_filtered_node_lenth = 0  # node这一串字符的长度
            for i in list_node:  # 形如 '1,li_36' '6,_22'  '7,_22'
                ori_filtered_node = i.split(',')[-1]  # li_36
                createVar['dict_'+i+'_allinfo'] = {}
                createVar['dict_'+i +
                          '_allinfo']['number'] = ori_filtered_node.split('_')[-1]
                if ori_filtered_node.startswith('_'):
                    # 孩子节点是叶子结点的父节点
                    createVar['dict_'+i+'_allinfo']['flag'] = 'child node'
                else:
                    createVar['dict_'+i +
                              '_allinfo']['flag'] = 'leaf node'  # 叶子节点
                '''下面这几行
                为了正确获取到
                字段一模一样的不同孩子节点（父节点）的 相关信息
                '''
                start_index = ori_filtered_node_index + \
                    ori_filtered_node_lenth  # 起始索引，第一次为0，每次更新到 剩余序列在整行中的绝对索引
                ori_index = ori_line[start_index:].index(
                    ori_filtered_node)  # 后半截序列里的node对应索引 后面进一步修改
                ori_filtered_node_index = ori_index+ori_filtered_node_index + \
                    ori_filtered_node_lenth  # 修改成在整行中的绝对索引
                createVar['dict_'+i +
                          '_allinfo']['index'] = ori_filtered_node_index
                createVar['dict_'+i +
                          '_allinfo']['depth'] = ori_line[:ori_filtered_node_index].count('(')-ori_line[:ori_filtered_node_index].count(')')+1
                # 对应深度层级=前半截序列左括号个数 - 右括号个数 + 1
                ori_filtered_node_lenth = len(ori_filtered_node)  # 节点长度
            # 判断基因家族是否收缩扩张
            for i in simplified_list_nodes_ids:  # li
                simplified_node_id = i
                i = list_node[simplified_list_nodes_ids.index(i)]  # 1,li_36
                current_number = createVar['dict_'+i+'_allinfo']['number']
                current_depth = createVar['dict_'+i+'_allinfo']['depth']
                '''下面几行的for循环
                为了找到当前节点深度-1的对应父节点
                以后有空再想想不用for循环的方式
                '''
                for j in list_node:  # 1,li_36
                    # 如果找到当前节点深度-1的对应父节点（child node是对整棵树而言，逻辑上是父节点）
                    if createVar['dict_'+j+'_allinfo']['depth'] == current_depth-1 and createVar['dict_'+j+'_allinfo']['flag'] == 'child node':
                        if current_number < createVar['dict_'+j+'_allinfo']['number']:
                            createVar['dict_'+simplified_node_id +
                                      '_allinfo']['shrunk'].append(family_name)
                        elif current_number == createVar['dict_'+j+'_allinfo']['number']:
                            createVar['dict_'+simplified_node_id +
                                      '_allinfo']['unchanged'].append(family_name)
                        elif current_number > createVar['dict_'+j+'_allinfo']['number']:
                            createVar['dict_'+simplified_node_id +
                                      '_allinfo']['expanded'].append(family_name)
    if trigger_flag == False:
        print('{}_rows_in_total_are_valid'.format(valid_row_count))

    print('\n---------------------------Step 03 check_the_number_of_gene_families------------------------------------')
    trigger_flag = False
    for simplified_node_id in simplified_list_nodes_ids:
        if valid_row_count != len(createVar['dict_'+simplified_node_id+'_allinfo']['shrunk']) +\
            len(createVar['dict_'+simplified_node_id + '_allinfo']['unchanged']) +\
                len(createVar['dict_'+simplified_node_id+'_allinfo']['expanded']):
            print("please_check_the_file")
            sys.exit(0)
    if trigger_flag == False:
        print('Done')

    print('\n---------------------------Step 04 write_file------------------------------------\n')
    with open(outtab_path, 'w') as out_handle:
        for simplified_node_id in simplified_list_nodes_ids:
            out_handle.write('node'.ljust(10)+'\tshrunk\n')
            out_handle.write(
                filtered_list_nodes_ids[simplified_list_nodes_ids.index(simplified_node_id)].ljust(10)+'\t')
            s = ''
            for i in createVar['dict_'+simplified_node_id+'_allinfo']['shrunk']:
                s += (i+'\t')
            out_handle.write(s.rstrip('\t')+'\n')
            out_handle.write('--'.ljust(10)+'\texpanded\n')
            out_handle.write('--'.ljust(10)+'\t')
            s = ''
            for i in createVar['dict_'+simplified_node_id+'_allinfo']['expanded']:
                s += (i+'\t')
            out_handle.write(s.rstrip('\t')+'\n\n')

    ###############################################################
    print('End Time : {}'.format(time.strftime(
        '%Y-%m-%d %H:%M:%S', time.localtime())))
    print('Already Run {}s'.format(time.time()-start_time))
    print('Done')
    ###############################################################
