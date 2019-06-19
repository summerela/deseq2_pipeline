import copy
from collections import OrderedDict

class Spreadsheet:
    def __init__(self,filename, header = False, delim='\t',string_columns = []):
        self.table = []
        self.delim = delim
        with open(filename) as f1:
            parsed_header = False
            for line in f1:
                if header and not parsed_header:
                    parsed_header = True
                    self.create_header(line)
                    continue
                data = line.strip().split(delim)
                if not data[0]:
                    continue
                for i,d in enumerate(data):
                    if i in string_columns:
                        continue
                    try:
                        fd = float(d)
                        data[i] = fd
                    except:
                        pass
                self.table.append(data)
        self.rows = len(self.table)
        if header:
            self.columns = max(len(self.table[0]),len(self.header.keys()))
        else:
            self.columns = len(self.table[0])
        for k in self.table:
            while len(k)<self.columns:
                k.append('')
    
    def create_header(self,line):
        headers = line.strip().split(self.delim)
        self.header = OrderedDict({})
        for i in range(len(headers)):
            self.header[headers[i]] = i
    
    def save(self,filename):
        output = []
        if hasattr(self,'header'):
            output.append(self.delim.join(self.header.keys()))
        for k in self.table:
            k = [str(i) for i in k]
            output.append(self.delim.join(k))
        with open(filename,'w') as f1:
            f1.write('\n'.join(output))
    
    def column(self,index):
        return [k[index] for k in self.table]
    
    def reduce_rows(self,keep):
        self.table = self.table[keep]
    
    def reduce_columns(self,keep):
        temp = copy.deepcopy(self.table)
        self.table =  [temp[k] for k in keep]
