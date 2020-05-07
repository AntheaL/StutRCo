import h5py
import reduce
from operator import iadd

def concatenate_dicts(dicts):
    res = dict()
    for key, value in dicts[0].items():
        if isinstance(value, dict):
            res[key] = dict()
            for k,v in value.items():
                res[key][k] = reduce(iadd, (d[key][k] for d in dicts))
        else:
            res[key] = reduce(iadd, (d[key] for d in dicts))

def save_h5(dst_path, data):
    with h5py.File(dst_path, "w") as dst:
        for key, value in info.items():
            if isinstance(value, dict):
                grp = dst.create_group(name=key)
                for k, v in value.items():
                    grp.create_dataset(name=k, data=v)
            else:
                dst.create_dataset(name=key, data=np.array(value))
