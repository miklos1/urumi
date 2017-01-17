from __future__ import absolute_import, print_function, division
from six import iteritems
from six.moves import range

import sys

_, filename = sys.argv

records = []

with open(filename) as f:
    for line in f.readlines():
        form, _1, dim, _2, degree, _3, nf, fc, time = line.split()
        records.append((form, int(dim), int(degree), int(nf), fc, float(time)))


def select(form, dim, nf, fc):
    for form_, dim_, degree, nf_, fc_, time in records:
        if form_ == form and dim == dim_ and nf == nf_ and fc == fc_:
            print("({}, {})".format(degree, time))


def plot(form, dim, fc, color):
    for nf in range(4):
        print(str.format(
            "\\addplot[color={color}!{opacity}!black,mark=*] coordinates {{",
            color=color,
            opacity=(100 - nf*20)
        ))
        select(form, dim, nf, fc)
        print("};")


colormap = {
    'uflacs': 'blue',
    # 'quadrature': 'red',
    'ffc-nonaffine': 'purple',
    'tsfc-quick': 'green',
    'tsfc-default': 'cyan',
    # 'tensor': 'orange',
}

for fc, color in iteritems(colormap):
    plot('helmholtz', 3, fc, color)
