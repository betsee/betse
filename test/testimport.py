#!/usr/bin/env python3

import yaml

data = ''

with open('testconfig.yaml', 'r') as yaml_file:
    data = yaml.load(yaml_file)


print(data)

a = 10.1e-10
b = 9.0e-10
c = 8.7e-10


data2 = {'Pmem_Na': a, 'Pmem_K': b, 'Pmem_Ca':c}


with open('testconfig.yaml', 'w') as yaml_file:
    yaml.dump(data2, yaml_file)



