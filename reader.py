import pysam


class Filestream:
    def __init__(self, file, start=None, n_iters=None):
        self.file = file
        self.len = None
        self.start = start
        self.n_iters = n_iters

    def __iter__(self):
        self.iter = 0
        with open(self.file) as f:
            for line in f:
                if self.start:
                    f.seek(self.start)
                line = self.line_function(line)
                if not line is None:
                    yield line
                self.iter += 1
                if self.n_iters and self.iter == self.n_iters:
                    break

    def __len__(self):
        if self.len is None:
            with open(self.file) as f:
                for i, l in enumerate(f):
                    pass
            self.len = i
        return self.len

    def line_function(self, line):
        return line


class FilestreamBED(Filestream):
    def line_function(self, _line):
        if _line.startswith("#"):
            return
        try:
            line = _line.strip().replace("chr", "").upper().split()
            region = {
                "chrom": line[0],#.split('_')[0],
                "start": int(line[1]),
                "stop": int(line[2]),
                "motif_len": int(line[3]),
            }
            if region["chrom"]=="00" or region["chrom"]=="0":
                raise ValueError(f"Invalid region {region} for line {line} stippred from {_line}")
            return region
        except (ValueError, IndexError):
            return


class FilestreamVCF(Filestream):
    def line_function(self, line):
        if line.startswith("#"):
            return
        try:
            line = line.strip().replace("chr", "").upper().split("\t")
            return {
                "chrom": line[0],
                "pos": int(line[1]),
                "id": line[2],
                "ref": line[3],
                "alt": line[4],
            }
        except (ValueError, IndexError):
            return
