import h5py
from functools import reduce
from operator import iadd
import numpy as np


def save_h5(dst_path, data):
    with h5py.File(dst_path, "w") as dst:
        for key, value in data.items():
            if isinstance(value, dict):
                grp = dst.create_group(name=key)
                for k, v in value.items():
                    grp.create_dataset(name=k, data=v)
            else:
                dst.create_dataset(name=key, data=np.array(value))


def load_h5(src_path):
    res = dict()
    with h5py.File(src_path, "r") as src:
        for key, value in src.items():
            if isinstance(value, h5py.Group):
                res[key] = dict()
                for k, v in value.items():
                    res[key][k] = np.array(v)
            else:
                res[key] = np.array(value)
    return res


def concatenate_files(dst_path, files):
    srcs = [h5py.File(f, "r") for f in files]
    with h5py.File(dst_path, "w") as dst:
        for key, value in srcs[0].items():
            if isinstance(value, h5py.Group):
                grp = dst.create_group(key)
                for k, v in value.items():
                    grp.create_dataset(
                        name=k,
                        data=np.concatenate([np.array(src[key][k]) for src in srcs]),
                    )
            else:
                dst.create_dataset(
                    name=key, data=np.concateate([d[key] for d in dicts])
                )
    [src.close() for src in srcs]


def concatenate_dicts(dicts):
    res = dict()
    for key, value in dicts[0].items():
        if isinstance(value, dict):
            res[key] = dict()
            for k, v in value.items():
                res[key][k] = np.concatenate([d[key][k] for d in dicts])
        else:
            res[key] = np.concatenate([src[key] for src in srcs])
    return res
