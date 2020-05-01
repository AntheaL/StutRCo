import pysam


class Filestream:
    def __init__(self, file):
        self.file = file
        self.len = None

    def __iter__(self):
        with open(self.file) as f:
            for line in f:
                line = self.line_function(line)
                if not line is None:
                    yield line

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
    def line_function(self, line):
        if line.startswith("#"):
            return
        try:
            line = line.strip().replace("chr", "").upper().split("\t")
            region = {
                "chrom": line[0],
                "start": int(line[1]),
                "stop": int(line[2]),
                "motif_len": int(line[3]),
                "num_ref_units": float(line[4]),
                "STR_id": line[5],
                "motif": line[6],
            }
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
